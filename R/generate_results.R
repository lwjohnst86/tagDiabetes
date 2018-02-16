#' (re-)generate the results of GEE, PLS, and correlation analyses.
#'
#' @return Creates several Rda files with lists of results in dataframes.
generate_results <- function() {
    generate_gee()
    generate_pls()
    generate_corr()
}

generate_gee <- function() {
    gee_unadj_df <- analyze_gee(covars = "YearsFromBaseline")
    gee_adj_df <- analyze_gee()
    gee_results <- list(
        unadj = gee_unadj_df,
        adj = gee_adj_df
    )
    devtools::use_data(gee_results, overwrite = TRUE)
}

generate_pls <- function() {
    # CV results are shown in pls.Rmd vignette
    pls_homa2_df <- analyze_pls(y = 'lHOMA2_S', cv = FALSE)
    pls_isi_df <- analyze_pls(y = 'lISI', cv = FALSE)
    pls_results <- list(
        homa2 = pls_homa2_df,
        isi = pls_isi_df
    )
    devtools::use_data(pls_results, overwrite = TRUE)
}

generate_corr <- function() {
    cor_conc_df <- analyze_corr(y = tg_conc)
    cor_pct_df <- analyze_corr(y = tg_pct)
    cor_tagfa_df <- analyze_corr_tagfa()
    cor_results <- list(
        conc = cor_conc_df,
        pct = cor_pct_df,
        tagfa = cor_tagfa_df
    )
    devtools::use_data(cor_results, overwrite = TRUE)
}
