
# Prepare data for QIC ----------------------------------------------------

#' Prepare the data for using in the QIC model comparisons.
#'
#' @param data Project data
#'
prep_qic_data <- function(data) {
    data %>%
        prep_gee_data() %>%
        dplyr::select(
            SID,
            VN,
            YearsFromBaseline,
            lISSI2,
            lISI,
            TotalTG,
            TotalNE,
            Sex,
            BiEthnicity,
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
}

# Clean QIC output --------------------------------------------------------

#' Clean data output from the QIC model selection functions.
#'
#' @param data QIC output data.
#'
#' @export
clean_qic <- function(data) {
    data %>%
        tibble::rownames_to_column() %>%
        dplyr::select(Model = rowname, QIC = IC, Delta = delta) %>%
        dplyr::mutate(id = seq_len(nrow(.)))
}

# Analyze using QIC -------------------------------------------------------

#' QIC output for insulin sensitivity.
#'
#' @param data Project data.
#'
#' @export
qic_is <- function(data = project_data) {
    # need to force into global environment because model.sel doesn't work
    qic_data <<- data %>%
        prep_qic_data() %>%
        dplyr::select(-lISSI2) %>%
        as.data.frame()

    M0 <- geepack::geeglm(
        lISI ~ TotalTG + YearsFromBaseline,
        data = qic_data,
        id = SID,
        family = stats::gaussian,
        corstr = 'ar1'
    )

    M1 <- stats::update(M0, . ~ . + TotalTG:YearsFromBaseline)
    M2 <- stats::update(M0, . ~ . + BaseAge + BiEthnicity + Sex)
    M3 <- stats::update(M0, . ~ . + BaseAge + BiEthnicity + Sex + Waist)
    M4 <- stats::update(M0, . ~ . + BaseAge + BiEthnicity + Sex + Waist + ALT)
    M5 <- stats::update(M0, . ~ . + BaseAge + BiEthnicity + Sex + Waist + ALT + MET)
    M6 <- stats::update(M0, . ~ . + BaseAge + BiEthnicity + Sex + Waist + ALT + MET + TotalNE)
    M7 <- stats::update(M0, . ~ . + BaseAge + BiEthnicity + Sex + Waist + ALT + MET + TotalNE +
                            AlcoholPerWk)
    M8 <- stats::update(M0, . ~ . + BaseAge + BiEthnicity + Sex + Waist + ALT + MET + TotalNE +
                            AlcoholPerWk + FamHistDiab)
    M9 <- stats::update(M0, . ~ . + BaseAge + BiEthnicity + Sex + Waist + ALT + MET + TotalNE +
                            AlcoholPerWk + FamHistDiab + TobaccoUse)

    MuMIn::model.sel(M0, M1, M2, M3, M4, M5, M6, M7, M8, M9,
                     rank = MuMIn::QIC) %>%
        clean_qic()
}

#' QIC output for beta-cell function.
#'
#' @param data Project data.
#'
#' @export
qic_bcf <- function(data = project_data) {
    # need to force into global environment because model.sel doesn't work
    qic_data <<- data %>%
        prep_qic_data() %>%
        dplyr::select(-lISI)

    M0 <- geepack::geeglm(
        lISSI2 ~ TotalTG + YearsFromBaseline,
        data = qic_data,
        id = SID,
        family = stats::gaussian,
        corstr = 'ar1'
    )

    M1 <- stats::update(M0, . ~ . + TotalTG:YearsFromBaseline)
    M2 <- stats::update(M0, . ~ . + BaseAge + BiEthnicity + Sex)
    M3 <- stats::update(M0, . ~ . + BaseAge + BiEthnicity + Sex + Waist)
    M4 <- stats::update(M0, . ~ . + BaseAge + BiEthnicity + Sex + Waist + ALT)
    M5 <- stats::update(M0, . ~ . + BaseAge + BiEthnicity + Sex + Waist + ALT + MET)
    M6 <- stats::update(M0, . ~ . + BaseAge + BiEthnicity + Sex + Waist + ALT + MET + TotalNE)
    M7 <- stats::update(M0, . ~ . + BaseAge + BiEthnicity + Sex + Waist + ALT + MET + TotalNE +
                            AlcoholPerWk)
    M8 <- stats::update(M0, . ~ . + BaseAge + BiEthnicity + Sex + Waist + ALT + MET + TotalNE +
                            AlcoholPerWk + FamHistDiab)
    M9 <- stats::update(M0, . ~ . + BaseAge + BiEthnicity + Sex + Waist + ALT + MET + TotalNE +
                            AlcoholPerWk + FamHistDiab + TobaccoUse)

    MuMIn::model.sel(M0, M1, M2, M3, M4, M5, M6, M7, M8, M9,
                     rank = MuMIn::QIC) %>%
        clean_qic()
}

#' Combine together the QIC datasets.
#'
#' @export
analyze_qic <- function() {
    qic_is_data <- qic_is() %>%
        dplyr::bind_rows(data.frame(Model = '**log(ISI)**', id = 0.5))

    qic_bcf_data <- qic_bcf() %>%
        dplyr::bind_rows(data.frame(Model = '**log(ISSI-2)**', id = 0.5)) %>%
        dplyr::mutate(id = max(id) + id)

    dplyr::bind_rows(qic_is_data, qic_bcf_data)
}

# Create table of QIC -----------------------------------------------------

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
