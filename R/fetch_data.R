#' Fetch data from the original source
#'
#' This function fetchs the main dataset, keeps variables relevant to
#' the analysis, restrict the sample size as needed, and lastly save
#' the new dataset as an `.RData` file.
#'
#' @return Saves the wrangled data into the data/ folder.
#' @export
#'
#' @examples
#' fetch_data()
#'
fetch_data <- function() {
    # Load the master dataset,
    ds.prep <- PROMISE::PROMISE_data %>%
        dplyr::filter(VN %in% c(1, 3, 6)) %>%
        ## Kick out Canoers
        dplyr::filter(is.na(Canoe)) %>%
        dplyr::tbl_df()

    print(paste0('Original dataset rows are ', dim(ds.prep)[1], ' and columns are ', dim(ds.prep)[2]))

    ##' Munge and wrangle the data into the final version.
    ds <- ds.prep %>%
        dplyr::select(
            SID, VN, BMI, Waist, HOMA, ISI, IGIIR, ISSI2, TAG, LDL, HDL, Chol,
            ALT, CRP, FamHistDiab, dplyr::matches('meds'), Age, Sex, Ethnicity,
            IFG, IGT, DM, MET, AlcoholPerWk, TobaccoUse, SelfEdu, Occupation,
            TotalTG, dplyr::matches('^tg\\d+'), Glucose0, Glucose120, TotalNE
        ) %>%
        dplyr::mutate(
            BaseTotalNE = TotalNE,
            BaseTotalTG = TotalTG,
            BaseTAG = ifelse(VN == 1, TAG, NA),
            lBaseTAG = log(BaseTAG),
            BaseAge = ifelse(VN == 1, Age, NA),
            FamHistDiab =
                plyr::mapvalues(FamHistDiab, c(0, 1:12),
                                c('No', rep('Yes', 12))) %>%
                as.factor(),
            AlcoholPerWk = plyr::mapvalues(
                as.factor(AlcoholPerWk),
                c('1', '2', '3', '4', '5', '6', '7'),
                c('1', '2', '3', '3', '3', '3', '3')
            ) %>%
                as.factor(),
            invHOMA = (1 / HOMA),
            linvHOMA = log(invHOMA),
            lISI = log(ISI),
            lIGIIR = log(IGIIR),
            lISSI2 = log(ISSI2),
            lALT = log(ALT),
            lTAG = log(TAG),
            MedsLipidsChol = ifelse(is.na(MedsLipidsChol), 0, MedsLipidsChol)
        ) %>%
        dplyr::arrange(SID, VN) %>%
        dplyr::group_by(SID) %>%
        tidyr::fill(TotalTG, dplyr::matches('^tg\\d+'), BaseTAG, lBaseTAG, BaseAge, TotalNE) %>%
        dplyr::ungroup()

    ds <- ds %>%
        dplyr::full_join(ds %>%
                      dplyr::filter(VN == 1) %>%
                      dplyr::mutate_each(dplyr::funs((. / TotalTG) * 100), dplyr::matches('^tg\\d+')) %>%
                      dplyr::select(SID, dplyr::matches('^tg\\d+')) %>%
                      stats::setNames(paste0('pct_', names(.))) %>%
                      dplyr::rename(SID = pct_SID),
                  by = 'SID')

    ds <- ds %>%
        dplyr::mutate(
            VN = plyr::mapvalues(VN, c(1, 3, 6), c(0, 1, 2)),
            f.VN = factor(VN, c(0, 1, 2), c('yr0', 'yr3', 'yr6')),
            Dysgly = plyr::mapvalues(as.character(IFG + IGT + DM), c('0', '1'), c('No', 'Yes')),
            Ethnicity =
                plyr::mapvalues(
                    Ethnicity,
                    c(
                        'African',
                        'European',
                        'First Nations',
                        'Latino/a',
                        'Other',
                        'South Asian'
                    ),
                    c('Other', 'European', 'Other', 'Latino/a',
                      'Other', 'South Asian')
                ),
            BiEthnicity = plyr::mapvalues(
                Ethnicity,
                c('Other', 'European', 'Latino/a', 'South Asian'),
                c('Non-European', 'European', 'Non-European', 'Non-European')
            )
        ) %>%
        dplyr::arrange(SID, VN) %>%
        dplyr::filter(!is.na(TotalTG))

    print(paste0('Working dataset rows are ', dim(ds)[1], ' and columns are ', dim(ds)[2]))

    # Final dataset object
    project_data <- ds

    # Save the dataset to the data/ folder.
    devtools::use_data(project_data, overwrite = TRUE)
    # Save the variable names as an internal dataset
    vars <- names(project_data)
    devtools::use_data(vars, internal = TRUE, overwrite = TRUE)
}
