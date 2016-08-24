# Plotting ----------------------------------------------------------------

#' Plots the distribution of the TAGFA composition.
#'
#' @param data Project data.
#'
#' @export
plot_tagfa <- function(data = project_data) {
    data %>%
        dplyr::filter(VN == 0) %>%
        dplyr::select(SID, dplyr::matches('^tg\\d+')) %>%
        tidyr::gather(Measure, Value,-SID) %>%
        stats::na.omit() %>%
        dplyr::mutate(Measure = renaming_fats(Measure)) %>%
        dplyr::mutate(
            order1 = substr(Measure, nchar(Measure), nchar(Measure)),
            order1 = ifelse(order1 == 0, 10, order1),
            order1 = ifelse(order1 == 'l', 20, order1),
            order1 = as.integer(order1)
        ) %>%
        dplyr::arrange(dplyr::desc(order1)) %>%
        dplyr::select(-order1) %>%
        dplyr::mutate(Measure = Measure %>%
                   factor(., levels = unique(.))) %>%
        seer::view_boxplots(
            'Measure', 'Value',
            xlab = 'Concentration (nmol/mL)',
            ylab = 'Triacylglyceride fatty acids') +
        graph_theme(ticks = FALSE)
}
