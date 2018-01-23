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
    ds.prep <- PROMISE.data::PROMISE %>%
        full_join(PROMISE.data::fattyacids %>%
                      select(SID, VN, matches("^tg|^TotalTG|^TotalNE"))) %>%
        filter(VN %in% c(0, 3, 6)) %>%
        ## Kick out Canoers
        filter(is.na(Canoe)) %>%
        tbl_df()

    print(paste0(
        'Original dataset rows are ',
        dim(ds.prep)[1],
        ' and columns are ',
        dim(ds.prep)[2]
    ))

    ##' Munge and wrangle the data into the final version.
    ds <- ds.prep %>%
        select(
            SID, VN, BMI, Waist, HOMA, HOMA2_S, ISI, IGIIR, ISSI2,
            TAG, LDL, HDL, Chol, YearsFromBaseline, Dysgly,
            ALT, CRP, FamHistDiab, matches('meds'), Age, Sex, Ethnicity,
            IFG, IGT, DM, MET, AlcoholPerWk, TobaccoUse, SelfEdu, Occupation,
            TotalTG, matches('^tg\\d+'), Glucose0, Glucose120, TotalNE
        ) %>%
        mutate(
            BaseTotalNE = TotalNE,
            BaseTotalTG = TotalTG,
            BaseTAG = ifelse(VN == 1, TAG, NA),
            BaseAge = ifelse(VN == 1, Age, NA),
            FamHistDiab =
                plyr::mapvalues(FamHistDiab, c(0, 1:12),
                                c('No', rep('Yes', 12))) %>%
                as.factor(),
            AlcoholPerWk = forcats::fct_other(
                AlcoholPerWk,
                keep = c("None", "<1 drink", "1-3 drinks"),
                other_level = "<3 drinks"),
            MedsLipidsChol = ifelse(is.na(MedsLipidsChol), 0, MedsLipidsChol)
        ) %>%
        arrange(SID, VN) %>%
        group_by(SID) %>%
        fill(TotalTG, matches('^tg\\d+'), BaseTAG, BaseAge, TotalNE) %>%
        ungroup()

    logged_vars <- ds %>%
        select(BaseTAG, HOMA2_S, ISI, IGIIR, ISSI2, ALT, TAG) %>%
        mutate_all(log) %>%
        rename_all(~paste0("l", .))

    pct_fats <- ds %>%
        filter(VN == 0) %>%
        select(SID, TotalTG, matches("^tg\\d+")) %>%
        mutate_at(vars(starts_with("tg")), funs((. / TotalTG) * 100)) %>%
        select(-TotalTG) %>%
        rename_at(vars(starts_with("tg")), funs(paste0("pct_", .)))

    ds <- ds %>%
        bind_cols(logged_vars) %>%
        full_join(pct_fats, by = 'SID')

    ds <- ds %>%
        mutate(
            VN = plyr::mapvalues(VN, c(0, 3, 6), c(0, 1, 2)),
            f.VN = factor(VN, c(0, 1, 2), c('yr0', 'yr3', 'yr6')),
            Dysgly = plyr::mapvalues(as.character(Dysgly), c("0", "1"), c("No", "Yes")),
            Ethnicity = forcats::fct_other(Ethnicity, keep = c("European", "Latino/a", "South Asian")),
            BiEthnicity = forcats::fct_other(Ethnicity, keep = "European", other_level = "Non-European")
        ) %>%
        arrange(SID, VN) %>%
        filter(!is.na(TotalTG))

    print(paste0('Working dataset rows are ', dim(ds)[1], ' and columns are ', dim(ds)[2]))

    # Final dataset object
    project_data <- ds

    # Save the dataset to the data/ folder.
    devtools::use_data(project_data, overwrite = TRUE)
    # Save the variable names as an internal dataset
    vars <- names(project_data)
    devtools::use_data(vars, internal = TRUE, overwrite = TRUE)
}
