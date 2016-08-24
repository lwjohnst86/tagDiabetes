# Functions for the PLS analysis
#
# Grab or combine data ----------------------------------------------------

#' Prepare the project data for PLS analysis.
#'
#' @param data Project data.
#' @param y Y or outcome variable.
#' @param x X or exposure variables of interest.
#'
#' @export
prep_pls_data <- function(data, y, x) {
    data %>%
        dplyr::filter(VN == 0) %>%
        dplyr::select_(.dots = c(y, x)) %>%
        stats::na.omit()
}

# Analyze -----------------------------------------------------------------

#' Compute PLS analysis.
#'
#' @param data Project data.
#' @param y Outcome variable.
#' @param x Exposure variables.
#' @param ncomp Number of components.
#'
#' @export
analyze_pls <- function(data = project_data,
                        y, x = tg_pct, ncomp = NULL) {
    data %>%
        prep_pls_data(y = y, x = x) %>%
        mason::design('pls') %>%
        mason::add_settings(ncomp = ncomp, validation = 'CV', cv.data = TRUE, cv.seed = 5436) %>%
        mason::add_variables('yvars', y) %>%
        mason::add_variables('xvars', x) %>%
        mason::construct() %>%
        mason::scrub()
}

# Plots -------------------------------------------------------------------

#' Plot the X loadings of the PLS results.
#'
#' @param data PLS output results.
#'
#' @export
plot_pls <- function(data) {
    seer::view_pls_xloadings(data, renaming.x = renaming_fats) +
        graph_theme(minor.grid.lines = FALSE)
}

# Calculations ------------------------------------------------------------

#' Calculate the correlation between the predicted and measured outcome.
#'
#' @param model The PLS output results.
#' @param test The test data to predict from.
#' @param ncomps The component to predict on.
#'
#' @export
calc_pred_corr <- function(model, test, ncomps = 1) {
    predicted <- predict(model, ncomp = ncomps, newdata = test)
    measured <- as.matrix(model.response(model.frame(formula(model),
                                                     data = test)))
    cor(predicted, measured)
}
