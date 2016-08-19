# Functions for the PLS analysis
#
# Grab or combine data ----------------------------------------------------

prep_pls_data <- function(data, y, x) {
    data %>%
        dplyr::filter(VN == 0) %>%
        dplyr::select_(.dots = c(y, x)) %>%
        stats::na.omit()
}

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
