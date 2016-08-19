# Functions for the GEE analysis
#
# Grab or combine data ----------------------------------------------------

#' Prepare the project data for analysis through GEE.
#'
#' @param data project data
#' @export
prep_gee_data <- function(data) {
    no_fattyacids <- data %>%
        dplyr::select(-dplyr::matches('pct_tg\\d+|^tg\\d+'),
                      -TotalNE,-TotalTG, -BaseTAG, -lBaseTAG)

    scaled_variables <- data %>%
        dplyr::filter(VN == 0) %>%
        dplyr::select(SID, TotalNE, TotalTG, BaseTAG, lBaseTAG,
                      dplyr::matches('pct_tg\\d+|^tg\\d+')) %>%
        dplyr::mutate_each(dplyr::funs(as.numeric(scale(.))), -SID)

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

# Analyze -----------------------------------------------------------------

#' Run GEE models on the prepared project data.
#'
#' @param data The project data
#' @param y outcomes (IS, BCF)
#' @param x predictors (TAGFA)
#' @param covariates to adjust for
#' @param intvar interaction variable
#' @param rename_x Function to rename x variables
#' @param rename_y Function to rename y variables
#' @export
analyze_gee <- function(data = project_data,
                        y = outcomes,
                        x = list(
                            tg_conc = tg_conc,
                            tg_pct = tg_pct,
                            tg_total = tg_totals
                        ),
                        covariates = c('VN', 'BaseAge', 'Sex', 'Ethnicity', 'Waist'),
                        intvar = NULL,
                        rename_x = renaming_fats,
                        rename_y = renaming_outcomes) {

    int <- !is.null(intvar)
    if (int) {
        extract_term <- ':'
    } else {
        extract_term <- 'Xterm$'
    }
    data %>%
        prep_gee_data() %>%
        mason::design('gee') %>%
        mason::add_settings(family = stats::gaussian(),
                            corstr = 'ar1', cluster.id = 'SID') %>%
        mason::add_variables('yvars', outcomes) %>%
        mason::add_variables('xvars', x[['tg_pct']]) %>%
        mason::add_variables('covariates', covariates) %>% {
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
        mason::construct() %>%
        mason::scrub() %>%
        mason::polish_filter(extract_term, 'term') %>%
        dplyr::mutate(unit = ifelse(grepl('pct', Xterms), 'mol%',
                                    ifelse(grepl('^tg\\d', Xterms), 'nmol/mL',
                                           'Totals'))) %>%
        mason::polish_transform_estimates(function(x) (exp(x) - 1) * 100) %>%
        mason::polish_renaming(rename_x, 'Xterms') %>%
        mason::polish_renaming(rename_y, 'Yterms') %>%
        dplyr::mutate(
            order1 = substr(Xterms, nchar(Xterms), nchar(Xterms)),
            order1 = ifelse(order1 == 0, 10, order1),
            order1 = ifelse(order1 == 'l', 20, order1),
            order1 = ifelse(order1 == 'G', 30, order1),
            order1 = as.integer(order1)
        ) %>%
        mason::polish_adjust_pvalue(method = 'BH') %>%
        dplyr::rename(unadj.p.value = p.value, p.value = adj.p.value) %>%
        dplyr::arrange(desc(order1)) %>%
        dplyr::mutate(Yterms = factor(
            Yterms,
            levels = c('log(1/HOMA-IR)', 'log(ISI)',
                       'log(IGI/IR)', 'log(ISSI-2)'),
            labels = c('log(1/HOMA-IR)', 'log(ISI)',
                       'log(IGI/IR)', 'log(ISSI-2)'),
            ordered = TRUE),
            Xterms = factor(Xterms, unique(Xterms))) %>%
        dplyr::select(-order1)

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
    gee_table_prep <- results %>%
        dplyr::mutate_each(dplyr::funs(format_rounding(., digits = digits)),
                           estimate, conf.low, conf.high) %>%
        dplyr::mutate(
            p.binary = ifelse(p.value <= 0.05, '\\*', ''),
            estimate.ci = paste0(estimate, ' (', conf.low, ', ', conf.high, ')', p.binary)
        ) %>%
        dplyr::select(Yterms, Xterms, unit, estimate.ci) %>%
        tidyr::spread(Yterms, estimate.ci)

    gee_table <-
        dplyr::bind_rows(
            data.frame(Xterms = paste0('**', unique(gee_table_prep$unit)[3], '**')),
            dplyr::filter(gee_table_prep, unit == 'Totals'),
            data.frame(Xterms = paste0('**', unique(gee_table_prep$unit)[1], '**')),
            dplyr::filter(gee_table_prep, unit == 'nmol/mL'),
            data.frame(Xterms = paste0('**', unique(gee_table_prep$unit)[2], '**')),
            dplyr::filter(gee_table_prep, unit == 'mol%')
        ) %>%
        dplyr::select(-unit) %>%
        dplyr::rename('Fatty acid' = Xterms)

    pander::pander(gee_table, missing = '', caption = caption)
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
            groups = 'unit~Yterms',
            legend.title = 'FDR-adjusted\np-value',
            xlab = 'Percent difference with 95% CI in the outcomes\nfor each SD increase in fatty acid',
            ylab = 'Triacylglycerol fatty acids'
            ) +
        graph_theme(ticks = FALSE)
}
