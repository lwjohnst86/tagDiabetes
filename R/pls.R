# Functions for the PLS analysis
#
# Grab, combine, or prep data ---------------------------------------------

#' Prepare the project data for PLS analysis.
#'
#' @param data Project data.
#' @param y Y or outcome variable.
#' @param x X or exposure variables of interest.
#'
#' @export
prep_pls_data <- function(data, y, x) {
    data %>%
        filter(VN == 0) %>%
        select_(.dots = c(y, x)) %>%
        na.omit()
}

#' Convert PLS model object to dataframe of loadings.
#'
#' @param .model PLS model object.
loadings_as_df <- function(.model) {
    loading_nums <- unclass(pls::loadings(.model))

    loading_nums %>%
        as_data_frame() %>%
        mutate(xvariables = dimnames(loading_nums)[[1]])

}

#' Add a column for "large" loadings (larger than abs(0.25)).
#'
#' @param .loadings PLS loadings dataset.
#'
append_large_loadings <- function(.loadings) {
    .loadings %>%
        mutate(xvariables = stringr::str_replace(xvariables, "pct_", "") %>%
                   PROMISE.misc::renaming_fa(keep.fraction = TRUE)) %>%
        gather(components, loadings, -xvariables) %>%
        mutate(loadings = as.numeric(loadings)) %>%
        group_by(components) %>%
        mutate(
            max_loading = max(loadings),
            min_loading = min(loadings),
            large_loadings = !between(loadings, -0.25, 0.25),
            xvariables = ifelse(large_loadings, xvariables, NA)
        ) %>%
        ungroup()
}

#' Extract scores from PLS results and include with original variables as a dataframe.
#'
#' @param model pls results
#'
scores_as_df <- function(model) {
    .scores <- pls::scores(model)
    attr(.scores, "explvar") <- NULL
    .scores <- .scores %>%
        unclass() %>%
        as.matrix() %>%
        as_tibble() %>%
        rename_all(funs(stringr::str_replace(., " ", "")))

    .data <- bind_cols(
        model$model$X %>%
            as_tibble(),
        model$model$Y %>%
            as_tibble()
        )

    bind_cols(.data, .scores)
}

# Analyze -----------------------------------------------------------------

#' Compute PLS analysis.
#'
#' @param data Project data.
#' @param y Outcome variable.
#' @param x Exposure variables.
#' @param ncomp Number of components.
#' @param cv Whether to use CV.
#'
#' @export
analyze_pls <- function(data = project_data,
                        y, x = tg_pct, ncomp = 4, cv = TRUE) {
    data %>%
        prep_pls_data(y = y, x = x) %>%
        design('pls') %>%
        add_settings(ncomp = ncomp, validation = 'CV', cv.data = cv, cv.seed = 5436) %>%
        add_variables('yvars', y) %>%
        add_variables('xvars', x) %>%
        construct() %>%
        scrub()
}

# Plots -------------------------------------------------------------------

#' Plots the loadings (weights/estimates) of the PLS model.
#'
#' @param .model PLS model object
#'
plot_pls_loadings <- function(.model) {
    total_expl_var <- pls::explvar(.model)[1:2] %>%
        sum() %>%
        round(1) %>%
        as.character()

    .model %>%
        loadings_as_df() %>%
        append_large_loadings() %>%
        filter(components %in% c("Comp 1", "Comp 2")) %>%
        mutate(components = stringr::str_replace(components, "omp ", "")) %>%
        ggplot(aes(x = loadings, y = components)) +
        geom_point(aes(alpha = large_loadings)) +
        ggrepel::geom_text_repel(
            aes(label = xvariables),
            size = 2.5,
            box.padding = 0.4,
            segment.alpha = 0.3,
            colour = "black"
        ) +
        labs(y = paste0("Components (", total_expl_var, "% total \nexplained variance)"),
             x = "Loading (coefficient) values for the component") +
        scale_alpha_discrete(guide = FALSE) +
        theme_classic()
}

#' Plots the PLS scores/predicted values from the model.
#'
#' @param .model PLS model object
#'
plot_pls_scores <- function(.model) {
    yvar <- dimnames(.model$model$Y)[[2]]
    .scores <- .model %>%
        scores_as_df()
    gather_vars <- grep("^Comp", names(.scores), value = TRUE)

    .output <- .scores %>%
        mutate_at(yvar, funs(as.factor(ntile(., 3)))) %>%
        select_at(vars(-starts_with("pct_"))) %>%
        gather_("Component", "Score", gather_vars) %>%
        mutate_at("Component", funs(stringr::str_replace(., "Comp", "Component ")))

    .output %>%
        ggplot(aes_string(x = yvar, y = "Score")) +
        geom_violin() +
        labs(y = "PLS scores for each individual and each component",
             x = "Tertiles of log(ISI)") +
        scale_x_discrete(labels = c("1st", "2nd", "3rd")) +
        facet_wrap( ~ Component)
}


#' Plot the X correlation loadings of the PLS results.
#'
#' @param data PLS output results.
#'
#' @export
plot_pls_corrcomps <- function(.model, outcome, title = NULL) {
    # xloadings <-
    expl_var_50 <- sqrt(1 / 2)

    fit <- cor(model.matrix(.model), pls::scores(.model)[, 1:2, drop = FALSE]) %>%
        as_tibble(rownames = "xvariables") %>%
        setNames(c('xvariables', 'C1', 'C2')) %>%
        mutate(xvariables = PROMISE.misc::renaming_fa(xvariables) %>%
                   stringr::str_replace("pct_", "")) %>%
        mutate(Fraction = factor(stringr::str_extract(xvariables, "NE|TG|PL|CE"))) %>%
        mutate(xvariables = ifelse(calc_radius(C1, C2) >= expl_var_50, xvariables, NA)) %>%
        mutate(large_loadings = ifelse(calc_radius(C1, C2) >= expl_var_50, TRUE, FALSE))

    circle_outer <- seer:::.circle_data(1)
    circle_inner <- seer:::.circle_data(sqrt(1 / 2))

    fig <- ggplot(fit, aes_string(x = "C1", y = "C2")) +
        geom_segment(aes(
            x = -1,
            y = 0,
            xend = 1,
            yend = 0
        ), colour = 'grey90') +
        geom_segment(aes(
            x = 0,
            y = -1,
            xend = 0,
            yend = 1
        ), colour = 'grey90') +
        geom_path(data = circle_outer, aes(x = x, y = y)) +
        geom_path(data = circle_inner, aes(x = x, y = y), linetype = 'dotted') +
        geom_point(data = fit, aes(alpha = large_loadings), show.legend = FALSE) +
        ggrepel::geom_text_repel(
            data = fit,
            aes(label = xvariables),
            size = 2.5,
            box.padding = 0.4,
            segment.alpha = 0.3
        ) +
        labs(
            x = paste0('C1 (', round(pls::explvar(.model)[1], 1), '% explained variance)'),
            y = paste0('C2 (', round(pls::explvar(.model)[2], 1), '% explained variance)'),
            title = title
        ) +
        theme_classic() +
        theme(legend.position = "none")
    fig
}


# Calculations ------------------------------------------------------------

#' Calculate the correlation between the predicted and actual outcome.
#'
#' @param model The PLS output results.
#' @param test The test data to predict from.
#' @param ncomps The component to predict on.
#'
#' @export
calc_pred_corr <- function(model, test, ncomps = 1) {
    predicted <- stats::predict(model, ncomp = ncomps, newdata = test)
    measured <- as.matrix(stats::model.response(
        stats::model.frame(formula(model), data = test))
        )
    corr <- broom::tidy(stats::cor.test(predicted, measured))[c(1, 3)]
    r <- aide::format_round(corr[1], 2)
    p <- aide::format_pval(corr[2])
    list(r = r, p = p, r_p = paste0('r=', r, ', p', p))

}

#' Calculates the distance from the center of two components.
#'
#' This is used to determine whether something is past the "threshold" of >50%
#' explained variance.
#'
#' @param x The x-axis variable, in this case Component 1
#' @param y The y-axis variable, in this case Component 2
calc_radius <- function(x, y) {
    sqrt(x ^ 2 + y ^ 2)
}

# large_contributors <- function(.model) {
#     expl_var_50 <- sqrt(1 / 2)
#     .model %>%
#         model.matrix() %>%
#         cor(pls::scores(.model)[, 1:2, drop = FALSE]) %>%
#         as_tibble(rownames = "fattyacid") %>%
#         rename_all(funs(stringr::str_replace(., " ", ""))) %>%
#         mutate(fattyacid = PROMISE.misc::renaming_fa(fattyacid) %>%
#                    stringr::str_replace("pct_", "")) %>%
#         filter(calc_radius(Comp1, Comp2) >= expl_var_50)
# }
