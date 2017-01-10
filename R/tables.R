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
        carpenter::add_rows(c('HOMA2_S', 'ISI', 'IGIIR', 'ISSI2'), carpenter::stat_medianIQR, digits = 1) %>%
        carpenter::add_rows(c('ALT', 'TAG', 'Chol', 'HDL', 'BaseTotalTG', 'BaseTotalNE', 'MET',
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

#' QIC table for manuscript.
#'
#' @param data Project data.
#' @param caption Caption for table.
#'
#' @export
table_qic_main <- function(data = project_data, caption = NULL) {
    prep_qic <- data %>%
        prep_gee_data() %>%
        dplyr::select(
            SID,
            VN,
            lISSI2,
            lISI,
            TotalTG,
            TotalNE,
            Sex,
            Ethnicity,
            MET,
            ALT,
            BaseAge,
            Waist,
            BMI,
            AlcoholPerWk,
            TobaccoUse,
            FamHistDiab
        ) %>%
        stats::na.omit()

    qic_is <- prep_qic %>%
        dplyr::select(-lISSI2)

    qic_bcf <- prep_qic %>%
        dplyr::select(-lISI)
    M0_ar_is <- geepack::geeglm(
        lISI ~ TotalTG + VN,
        data = qic_is,
        id = SID,
        family = stats::gaussian,
        corstr = 'ar1'
    )
    M0 <- M0_ar_is
    Full <- stats::update(M0, . ~ . + Waist + BaseAge + Ethnicity + Sex + ALT + MET +
                AlcoholPerWk + FamHistDiab + TobaccoUse + TotalNE)
    M1 <- stats::update(M0, . ~ . + Waist + BaseAge + Ethnicity + Sex)
    M2 <- stats::update(M0, . ~ . + Waist + BaseAge + Ethnicity + Sex + ALT)
    M3 <- stats::update(M0, . ~ . + Waist + BaseAge + Ethnicity + Sex + ALT + MET)
    M4 <- stats::update(M0, . ~ . + Waist + BaseAge + Ethnicity + Sex + ALT + MET + AlcoholPerWk)
    M5 <- stats::update(M0, . ~ . + Waist + BaseAge + Ethnicity + Sex + ALT + MET + AlcoholPerWk +
                          FamHistDiab)
    M6 <- stats::update(M0, . ~ . + Waist + BaseAge + Ethnicity + Sex + ALT + FamHistDiab + TobaccoUse)
    M7 <- stats::update(M0, . ~ . + Waist + BaseAge + Ethnicity + Sex + ALT + FamHistDiab)
    M8 <- stats::update(M0, . ~ . + Waist + BaseAge + Ethnicity + Sex + ALT + MET + FamHistDiab)
    M9 <- stats::update(M0, . ~ . + Waist + BaseAge + Ethnicity + Sex + ALT + MET + TotalNE)
    M10 <- stats::update(M0, . ~ . + Waist + BaseAge + Ethnicity + Sex + ALT + MET + TotalNE +
                          FamHistDiab)

    mod_is <-
        MuMIn::model.sel(M0, Full, M1, M2, M3, M4, M5, M6, M7, M8, M9, M10,
                         rank = MuMIn::QIC) %>%
        clean_qic() %>%
        dplyr::bind_rows(data.frame(Model = '**log(ISI)**', id = 0.5))

    M0_ar_bcf <- geepack::geeglm(
        lISSI2 ~ TotalTG + VN,
        data = qic_bcf,
        id = SID,
        family = stats::gaussian,
        corstr = 'ar1'
    )
    M0 <- M0_ar_bcf
    Full <- stats::update(M0, . ~ . + Waist + BaseAge + Ethnicity + Sex + ALT + MET +
                AlcoholPerWk + FamHistDiab + TobaccoUse + TotalNE)
    M1 <- stats::update(M0, . ~ . + Waist + BaseAge + Ethnicity + Sex)
    M2 <- stats::update(M0, . ~ . + Waist + BaseAge + Ethnicity + Sex + ALT)
    M3 <- stats::update(M0, . ~ . + Waist + BaseAge + Ethnicity + Sex + ALT + MET)
    M4 <- stats::update(M0, . ~ . + Waist + BaseAge + Ethnicity + Sex + ALT + MET + AlcoholPerWk)
    M5 <- stats::update(M0, . ~ . + Waist + BaseAge + Ethnicity + Sex + ALT + MET + AlcoholPerWk +
                          FamHistDiab)
    M6 <- stats::update(M0, . ~ . + Waist + BaseAge + Ethnicity + Sex + ALT + FamHistDiab + TobaccoUse)
    M7 <- stats::update(M0, . ~ . + Waist + BaseAge + Ethnicity + Sex + ALT + FamHistDiab)
    M8 <- stats::update(M0, . ~ . + Waist + BaseAge + Ethnicity + Sex + ALT + MET + FamHistDiab)
    M9 <- stats::update(M0, . ~ . + Waist + BaseAge + Ethnicity + Sex + ALT + MET + TotalNE)
    M10 <- stats::update(M0, . ~ . + Waist + BaseAge + Ethnicity + Sex + ALT + MET + TotalNE +
                          FamHistDiab)

    mod_bcf <-
        MuMIn::model.sel(M0, Full, M1, M2, M3, M4, M5, M6, M7, M8, M9, M10,
                         rank = MuMIn::QIC) %>%
        clean_qic() %>%
        dplyr::mutate(id = max(id) + id) %>%
        {
            dplyr::bind_rows(., data.frame(Model = '**log(ISSI-2)**', id = min(.['id']) - 0.5))
        }

    dplyr::bind_rows(mod_is, mod_bcf) %>%
        table_qic(caption = caption)

}
