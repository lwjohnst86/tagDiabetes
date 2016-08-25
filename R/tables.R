# Tables ------------------------------------------------------------------

#' Create a table of the basic characteristics of PROMISE.
#'
#' @param data Project data.
#' @param caption Table caption.
#'
#' @export
table_basic <- function(data = project_data, caption = NULL) {
    data %>%
        dplyr::mutate(Ethnicity = ifelse(VN == 0, as.character(Ethnicity), NA),
                      Ethnicity = as.factor(Ethnicity),
                      Sex = ifelse(VN == 0, as.character(Sex), NA),
                      Sex = as.factor(Sex)) %>%
        carpenter::outline_table('f.VN') %>%
        carpenter::add_rows(c('HOMA', 'ISI', 'IGIIR', 'ISSI2'), carpenter::stat_medianIQR, digits = 1) %>%
        carpenter::add_rows(c('ALT', 'TAG', 'Chol', 'HDL', 'BaseTotalNE', 'MET',
                              'Age', 'BMI', 'Waist'),
                            carpenter::stat_meanSD, digits = 1) %>%
        carpenter::add_rows(c('Ethnicity', 'Sex'), carpenter::stat_nPct, digits = 0) %>%
        carpenter::renaming('rows', renaming_table_rows) %>%
        carpenter::renaming("header", c('Measure', 'Baseline', '3-yr', '6-yr')) %>%
        carpenter::build_table(caption = caption)
}

#' Create a table of the values for the conc and mol% of TAGFA.
#'
#' @param data Project data.
#' @param caption Table caption.
#'
#' @export
table_tagfa <- function(data = project_data, caption = NULL) {
    fa <- function(pattern, x = c(tg_conc, 'TotalTG'))
        grep(pattern, x, value = TRUE)
    tgfa <- c(fa('3$'), fa('6$'), fa('7$'), fa('9$'), fa('0$'), fa('TotalTG$'))

    data %>%
        dplyr::filter(VN == 0) %>%
        tidyr::gather(FA, Value, dplyr::matches('tg\\d')) %>%
        dplyr::mutate(unit = as.factor(ifelse(grepl('pct_', FA), 'mol%', 'nmol/mL')),
                      FA = gsub('pct_', '', FA)) %>%
        tidyr::spread(FA, Value) %>%
        dplyr::mutate(TotalTG = ifelse(unit == 'mol%', NA, TotalTG),
                      unit = factor(unit, levels = c('nmol/mL', 'mol%'), ordered = TRUE)) %>%
        carpenter::outline_table('unit') %>%
        carpenter::add_rows(tgfa, carpenter::stat_meanSD, digits = 1) %>%
        carpenter::renaming("header", c('TAGFA', 'Concentrations (nmol/mL)', 'Proportion (mol%)')) %>%
        carpenter::renaming('rows', renaming_fats) %>%
        carpenter::build_table(caption = caption)
}

table_corr <- function(results, caption = NULL) {
    results %>%
        tidyr::spread(Vars2, Correlations) %>%
        dplyr::rename(TAGFA = Vars1) %>%
        pander::pander()
}

#' Create a table of the QIC from the QIC cleaned output.
#'
#' @param data Cleaned QIC output.
#' @param caption Table caption.
#'
#' @export
table_qic <- function(data, caption = NULL) {
    data %>%
        dplyr::arrange(id) %>%
        dplyr::select(-id) %>%
        dplyr::mutate(QIC = format_rounding(QIC, 1),
                      Delta = format_rounding(Delta, 1)) %>%
        pander::pander(missing = '', caption = caption)
}
