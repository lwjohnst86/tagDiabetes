# Analyze -----------------------------------------------------------------

analyze_gee <- function(data = project_data,
                        y = outcomes,
                        x = list(
                            tg_conc = tg_conc,
                            tg_pct = tg_pct,
                            tg_total = tg_total
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

analyze_pls <- function(data = project_data,
                        y, x = tg_pct, ncomp = 2, type = c('train', 'test')) {
    data %>%
        prep_pls_data(y = y, x = x) %>%
        mason::add_settings(ncomp = ncomp, cv.data = TRUE, cv.seed = 43145) %>%
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
