# Functions to run the correlation analysis
#
# Analyze -----------------------------------------------------------------

#' Correlation analysis.
#'
#' @param data Project data
#' @param x The covariates and outcomes
#' @param y The fatty acids.
#'
#' @export
analyze_corr <-
    function(data = project_data,
             x = c(outcomes, 'BMI', 'Waist', 'Age', 'lALT', 'lTAG', 'Chol', 'HDL', 'LDL'),
             y = tg_conc) {

    data %>%
        dplyr::filter(VN == 0) %>%
        mason::design('cor') %>%
        mason::add_settings(method = 'pearson', use = 'complete.obs') %>%
        mason::add_variables('yvars', y) %>%
        mason::add_variables('xvars', x) %>%
        mason::construct() %>%
        mason::scrub() %>%
        mason::polish_renaming(renaming_fats, 'Vars2') %>%
        mason::polish_renaming(function(x)
            gsub('l(ALT|TAG|IGIIR|invHOMA|ISI|ISSI2)', '\\1', x) %>%
                renaming_outcomes(), 'Vars1') %>%
        dplyr::mutate(
            order1 = substr(Vars2, nchar(Vars2), nchar(Vars2)),
            order1 = ifelse(order1 == 0, 10, order1),
            order1 = ifelse(order1 == 'l', 20, order1),
            order1 = as.integer(order1)
        ) %>%
        dplyr::arrange(dplyr::desc(order1)) %>%
        dplyr::select(-order1) %>%
        dplyr::mutate(Vars2 = factor(Vars2, unique(Vars2)),
                      Vars1 = factor(Vars1, unique(Vars1)),
                      Correlations = round(Correlations, 2))
    }

# Plotting ----------------------------------------------------------------

#' Correlation heatmap plot.
#'
#' @param results Correlation results
#'
#' @export
plot_heatmap <- function(results) {
     results %>%
        seer::view_heatmap(
            y = 'Vars2',
            x = 'Vars1',
            ylab = 'Triacylglycerol fatty acids (nmol/mL)',
            number.colours = 5,
            values.size = 4) +
        graph_theme(ticks = FALSE, legend.pos = 'right')
}
