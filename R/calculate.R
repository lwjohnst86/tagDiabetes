# Calculate or extract for inline -----------------------------------------

#' Calculate the number of individuals who converted to pre-DM or DM.
#'
#' @param data Data with conversion numbers.
#'
#' @export
calc_convert_dys <- function(data) {
    conversion_dm <-
        table(data['ConvertDM'])[2] %>%
        {
            paste0(., ' (', round((. / 477) * 100, 0), '%)')
        }

    conversion_pre <-
        table(data['ConvertPreDM'])[2] %>%
        {
            paste0(., ' (', round((. / 477) * 100, 0), '%)')
        }

    list(pre = conversion_pre, dm = conversion_dm)
}

#' Calculate the percent contribution of each TAGFA of the whole fraction.
#'
#' @param data Project data.
#'
#' @export
calc_tagfa_percent <- function(data = project_data) {
    over_10pct_tagfa <- data %>%
        dplyr::filter(VN == 0) %>%
        dplyr::select(dplyr::matches('pct_tg')) %>%
        tidyr::gather(fat, value) %>%
        dplyr::group_by(fat) %>%
        dplyr::summarise(pct = mean(value, na.rm = TRUE)) %>%
        dplyr::mutate(fat = renaming_fats(fat),
               c = paste0(fat, ' (', aide::format_round(pct, 1), '%)')) %>%
        dplyr::arrange(dplyr::desc(pct)) %>%
        dplyr::filter(pct >= 10)

    return(over_10pct_tagfa)
}

#' Calculate the percent change over time for the outcome variables.
#'
#' @param data Project data.
#'
#' @export
calc_outcome_changes <- function(data = project_data) {
    prep.data <- data %>%
        dplyr::select(f.VN, HOMA2_S, ISI, IGIIR, ISSI2) %>%
        tidyr::gather(Measure, Value,-f.VN) %>%
        stats::na.omit() %>%
        dplyr::group_by(Measure, f.VN) %>%
        dplyr::summarise(med = median(Value),
                         n = n()) %>%
        dplyr::ungroup()

    sample_size <- prep.data$n %>%
        {paste0(min(.), '-', max(.))}

    change_over_time <- prep.data %>%
        dplyr::select(-n) %>%
        tidyr::spread(f.VN, med) %>%
        dplyr::mutate(pctChg = ((yr6 - yr0) / yr0) * 100) %>%
        dplyr::select(pctChg) %>%
        abs() %>%
        round(0) %>%
        {paste0(min(.), '% to ', max(.), '%')}

    pval <- mason::design(data, 'gee') %>%
        mason::add_settings(family = stats::gaussian(), corstr = 'ar1', cluster.id = 'SID') %>%
        mason::add_variables('yvars', c('lHOMA2_S', 'lISI', 'lIGIIR', 'lISSI2')) %>%
        mason::add_variables('xvars', 'VN') %>%
        mason::construct() %>%
        mason::scrub() %>%
        mason::polish_filter('Xterm$', 'term') %>%
        dplyr::summarise(p.value = mean(p.value)) %>%
        dplyr::mutate(p.value = aide::format_pval(p.value))

    change_outcomes <- list(n = sample_size, chg = change_over_time, p = pval)

    return(change_outcomes)
}


#' Calculate the n and % of total for a given discrete variable (e.g. Sex).
#'
#' @param x The discrete variable data.
#' @param group Which group to calculate n and % for.
#'
#' @export
calc_discr_npct <- function(x, group = c('Yes', 'Female', 'European')) {
    nums <- table(x)
    pct <- (nums[group]/sum(nums)) * 100
    data.frame(
        n = nums[group],
        pct = paste0(round(pct, 1), '%'),
        npct = paste0(nums[group], ' (', paste0(round(pct, 1), '%)'))
    )
}

calc_gee_ave <- function(results) {
    results %>%
        dplyr::filter(p.value <= 0.05) %>%
        dplyr::group_by(Yterms) %>%
        dplyr::summarise(b = aide::format_round(abs(mean(estimate))))
}


#' Extract the estimates and confidence interval for the GEE results.
#'
#' @param results GEE analysis results
#'
#' @export
calc_gee_estci <- function(results) {
    b_ci <- results %>%
        dplyr::mutate(
            est = aide::format_round(estimate),
            ci = paste0(
                aide::format_round(conf.low),
                ' to ',
                aide::format_round(conf.high)
            ),
            eci = paste0('(beta: ', est, ', CI: ', ci, ')')
        ) %>%
        dplyr::filter(p.value <= 0.05) %>%
        dplyr::arrange(Yterms, dplyr::desc(estimate)) %>%
        dplyr::select(Yterms, Xterms, est, eci)

    outcome_b <- results %>%
        dplyr::filter(p.value <= 0.05) %>%
        dplyr::mutate(dir = ifelse(estimate < 0, 'Neg', 'Pos')) %>%
        dplyr::group_by(Yterms, unit, dir) %>%
        dplyr::summarise(rng = aide::min_max(estimate))

    list(b_ci = b_ci, outcome_b = outcome_b)
}

#' Extracts and computes various information from the correlation results.
#'
#' @param data Correlation results data.
#'
#' @export
calc_cor <- function(data = cor_df) {
    data %>%
        {
            high_cor <- dplyr::filter(., !dplyr::between(Correlations, -0.3, 0.3)) %>%
                dplyr::mutate(Direction = dplyr::if_else(Correlations > 0.3, 'Positive', 'Negative'))

            list(
                rng = high_cor %>%
                    dplyr::group_by(Vars1, Direction) %>%
                    dplyr::summarize(rng = paste0(min(Correlations), ' to ', max(Correlations))) %>%
                    dplyr::ungroup() %>%
                    tidyr::spread(Vars1, rng),
                pos_fa = dplyr::filter(high_cor, Direction == 'Positive') %>%
                    dplyr::arrange(dplyr::desc(Correlations)),
                neg_fa = dplyr::filter(high_cor, Direction == 'Negative') %>%
                    dplyr::arrange(Correlations)
            )
        }
}

#' Calculate the sample size at each visit collection.
#'
#' May not be needed.
#'
#' @param data Project data
#'
#' @export
calc_n_at_each_visit <- function(data = project_data) {
    data %>%
        dplyr::select(SID, VN) %>%
        stats::na.omit() %>%
        dplyr::group_by(VN) %>%
        dplyr::summarise(n = n())
}

#' Calculate the number of participants who attended one, two, or three visits.
#'
#' @param data Project data.
#'
#' @export
calc_n_for_visits <- function(data = project_data) {
    data %>%
        dplyr::select(SID, VN, f.VN) %>%
        stats::na.omit() %>%
        tidyr::spread(f.VN, VN) %>%
        dplyr::mutate(Visits = paste(yr0, yr3, yr6, sep = '-')) %>%
        .$Visits %>%
        table()
}

#' Mean and SD followup time for those who attended the 6 year visit.
#'
#' @param data Project data
#'
#' @export
calc_followup_time <- function(data = project_data) {
    data %>%
        dplyr::arrange(SID, VN) %>%
        dplyr::group_by(SID) %>%
        dplyr::slice(n()) %>%
        dplyr::ungroup() %>%
        dplyr::summarise(MeanFollowup = aide::ave_sd(YearsFromBaseline))
}

#' Percent of participants who attended each of the 3 visits.
#'
#' @export
calc_pct_full_visits <- function() {
    pct_val <- calc_n_for_visits()[1] / sum(calc_n_for_visits())
    aide::format_round(pct_val * 100)
}
