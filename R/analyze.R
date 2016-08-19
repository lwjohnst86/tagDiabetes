# Analyze -----------------------------------------------------------------

analyze_pls <- function(data = project_data,
                        y, x = tg_pct, ncomp = 2,
                        type = c('train', 'test', 'full')) {
    data %>%
        prep_pls_data(y = y, x = x) %>%
        mason::design('pls') %>%
        mason::add_settings(ncomp = ncomp, cv.data = TRUE, cv.seed = 12345) %>%
        mason::add_variables('yvars', y) %>%
        mason::add_variables('xvars', x) %>%
        mason::construct(cv = type) %>%
        mason::scrub()
}

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
        mason::polish_renaming(renaming_fa, 'Vars2') %>%
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
