# Functions for the GEE analysis
#
# Grab or combine data ----------------------------------------------------

#' Prepare the project data for analysis through GEE.
#'
#' @param data project data
#' @export
prep_gee_data <- function(data) {
    data <- data %>%
        dplyr::filter(DM != 1)

    no_fattyacids <- data %>%
        dplyr::select(-dplyr::matches('pct_tg\\d+|^tg\\d+'),
                      -TotalNE,-TotalTG, -BaseTAG, -lBaseTAG)

    scaled_variables <- data %>%
        dplyr::filter(VN == 0) %>%
        dplyr::select(SID, TotalNE, TotalTG, BaseTAG, lBaseTAG,
                      dplyr::matches('pct_tg\\d+|^tg\\d+')) %>%
        dplyr::mutate_at(dplyr::vars(-SID), dplyr::funs(as.numeric(scale(.))))

    dplyr::full_join(
            no_fattyacids,
            scaled_variables,
            by = 'SID'
        ) %>%
        dplyr::group_by(VN) %>%
        dplyr::mutate(
            Waist = as.numeric(scale(Waist)),
            ALT = as.numeric(scale(ALT)),
            BaseAge = as.numeric(scale(BaseAge)),
            MET = as.numeric(scale(MET))
        ) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(SID, VN)
}

#' Add the extracted PLS scores to the data for use in GEE.
#'
#' @param data The project data.
#' @param pls_results Results dataset from PLS.
#' @export
prep_gee_pls_data <- function(data = project_data, pls_results) {
    sid_vn <- data %>%
        dplyr::filter(VN == 0) %>%
        dplyr::select_(.dots = c(outcomes, tg_pct, 'SID', covariates)) %>%
        stats::na.omit() %>%
        dplyr::full_join(dplyr::select(project_data, SID, VN)) %>%
        dplyr::select_(.dots = c(outcomes, 'SID', 'VN'))

    pls_results %>%
    {
        dplyr::bind_cols(tibble::as_data_frame(.$model$Y),
                         tibble::as_data_frame(unclass(.$scores))[1:2])
    } %>%
        dplyr::rename(Comp1 = `Comp 1`, Comp2 = `Comp 2`) %>%
        dplyr::full_join(sid_vn) %>%
        dplyr::select(SID, VN, Comp1, Comp2) %>%
        dplyr::full_join(data) %>%
        dplyr::arrange(SID, VN) %>%
        tidyr::fill(Comp1, Comp2)
}

# Analyze -----------------------------------------------------------------

#' Run GEE models on the prepared project data.
#'
#' @param data The project data
#' @param y outcomes (IS, BCF)
#' @param x predictors (TAGFA)
#' @param intvar interaction variable
#' @param rename_x Function to rename x variables
#' @param covars Covariates to include in the model
#' @param data_prep Whether to clean the dataset before hand (either basic prep
#'   or extracting PLS scores).
#' @param pls_results Optional; if doing the data_prep for PLS, need to include
#'   the pls results dataset.
#' @param rename_y Function to rename y variables
#'
#' @export
analyze_gee <- function(data = project_data,
                        y = outcomes,
                        x = list(
                            tg_conc = tg_conc,
                            tg_pct = tg_pct,
                            tg_total = tg_totals
                        ),
                        covars = covariates,
                        intvar = NULL,
                        rename_x = renaming_fats,
                        rename_y = renaming_outcomes,
                        data_prep = c('basic', 'pls', 'none'),
                        pls_results = NULL) {

    int <- !is.null(intvar)
    if (int) {
        extract_term <- ':'
    } else {
        extract_term <- 'Xterm$'
    }

    data_prep <- match.arg(data_prep)
    switch(data_prep,
           none = {data <- data},
           basic = {data <- prep_gee_data(data)},
           pls = {data <- prep_gee_pls_data(data, pls_results = pls_results)})

    data %>%
        mason::design('gee') %>%
        mason::add_settings(family = stats::gaussian(),
                            corstr = 'ar1', cluster.id = 'SID') %>%
        mason::add_variables('yvars', y) %>%
        {if (is.list(x)) {
            mason::add_variables(., 'xvars', x[['tg_pct']]) %>%
                mason::add_variables('covariates', covars) %>% {
                    if (int) {
                        mason::add_variables(., 'interaction', intvar)
                    } else {
                        .
                    }
                } %>%
                mason::construct() %>%
                mason::add_variables('xvars', x[['tg_conc']]) %>%
                mason::construct() %>%
                mason::add_variables('xvars', x[['tg_total']]) %>%
                mason::construct()
        } else {
            mason::add_variables(., 'xvars', x) %>%
                mason::construct()
        }} %>%
        mason::scrub() %>%
        mason::polish_filter(extract_term, 'term') %>%
        categorize_fa_model_output() %>%
        transform_continuous_model_output() %>%
        rename_terms_model_output(rename_x = rename_x, rename_y = rename_y) %>%
        padjust_model_output() %>%
        order_model_output() %>%
        set_outcome_order_model_output()

}

categorize_fa_model_output <- function(gee_results) {
    gee_results %>%
        dplyr::mutate(unit = ifelse(
            grepl('pct', Xterms),
            'mol%',
            ifelse(grepl('^tg\\d', Xterms), 'nmol/mL',
                   'Totals')
        ))
}

transform_continuous_model_output <- function(gee_results) {
    gee_results %>%
        mason::polish_transform_estimates(function(x) (exp(x) - 1) * 100)
}

rename_terms_model_output <- function(gee_results, rename_x, rename_y) {
    gee_results %>%
        dplyr::mutate_at('Xterms', dplyr::funs(rename_x)) %>%
        dplyr::mutate_at('Yterms', dplyr::funs(rename_y))
}

padjust_model_output <- function(gee_results) {
    gee_results %>%
        mason::polish_adjust_pvalue(method = 'BH') %>%
        dplyr::rename(unadj.p.value = p.value, p.value = adj.p.value)
}

order_model_output <- function(gee_results) {
    gee_results %>%
        dplyr::mutate(
            order1 = substr(Xterms, nchar(Xterms), nchar(Xterms)),
            order1 = ifelse(order1 == 0, 10, order1),
            order1 = ifelse(order1 == 'l', 20, order1),
            order1 = ifelse(order1 == 'G', 30, order1),
            order1 = as.integer(order1)
        ) %>%
        dplyr::arrange(desc(order1)) %>%
        dplyr::mutate(Xterms = factor(Xterms, unique(Xterms))) %>%
        dplyr::select(-order1)
}

set_outcome_order_model_output <- function(gee_results) {
    gee_results %>%
        dplyr::mutate(Yterms = factor(
            Yterms,
            levels = c('log(HOMA2-IR)', 'log(HOMA2-%S)', 'log(ISI)',
                       'log(IGI/IR)', 'log(ISSI-2)'),
            labels = c('log(HOMA2-IR)', 'log(HOMA2-%S)', 'log(ISI)',
                       'log(IGI/IR)', 'log(ISSI-2)'),
            ordered = TRUE))
}

# Tables ------------------------------------------------------------------

#' Table to present the main GEE analysis.
#'
#' @param results Result data frame from the GEE analysis
#' @param caption Caption for the table
#' @param digits Digits to show for the table
#'
#' @export
table_gee_main <- function(results, caption = NULL, digits = 1) {
    table_data <- results %>%
        dplyr::mutate_at(
            dplyr::vars(estimate, conf.low, conf.high),
            dplyr::funs(format_rounding(., digits = digits))
        ) %>%
        dplyr::mutate(
            p.binary = ifelse(p.value <= 0.05, '\\*', ''),
            estimate.ci = paste0(estimate, ' (', conf.low, ', ', conf.high, ')',
                                 p.binary)
        ) %>%
        dplyr::select(Yterms, Xterms, unit, estimate.ci) %>%
        tidyr::spread(Yterms, estimate.ci) %>%
        dplyr::mutate(Xterms = as.character(Xterms))

    dplyr::bind_rows(
            tibble::data_frame(Xterms = paste0('**Totals**')),
            dplyr::filter(table_data, unit == 'Totals'),
            tibble::data_frame(Xterms = paste0('**nmol/mL**')),
            dplyr::filter(table_data, unit == 'nmol/mL'),
            tibble::data_frame(Xterms = paste0('**mol%**')),
            dplyr::filter(table_data, unit == 'mol%')
        ) %>%
        dplyr::select(-unit) %>%
        dplyr::rename('Fatty acid' = Xterms) %>%
        as.data.frame() %>%
        pander::pander(missing = '', caption = caption)
}

# Plotting ----------------------------------------------------------------

#' Plot the GEE results in a Forest plot style.
#'
#' @param results Results data frame from the GEE analysis
#'
#' @export
plot_gee_main <- function(results) {
    results %>%
        seer::view_main_effect(
            graph.options = 'dot.size',
            dot.colour = "black",
            groups = 'unit~Yterms',
            legend.title = 'FDR-adjusted\np-value',
            xlab = 'Percent difference with 95% CI in the outcomes\nfor each SD increase in fatty acid',
            ylab = 'Triacylglycerol fatty acids'
            ) +
        graph_theme(ticks = FALSE) +
        ggplot2::facet_grid(unit~Yterms, switch = 'y',
                            scales = 'free_y',
                            space = 'free_y')
}
