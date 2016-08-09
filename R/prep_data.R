# Grab or combine data ----------------------------------------------------

prep_pls_data <- function(data, y, x) {
    data %>%
        dplyr::filter(VN == 0) %>%
        dplyr::select_(y, x) %>%
        stats::na.omit()
}

prep_gee_data <- function(data) {
    no_fattyacids <- data %>%
        dplyr::select(-dplyr::matches('pct_tg\\d+|^tg\\d+'),
                      -TotalNE,-TotalTG, -BaseTAG, -lBaseTAG)

    scaled_variables <- data %>%
        dplyr::filter(VN == 0) %>%
        dplyr::select(SID, TotalNE, TotalTG, BaseTAG, lBaseTAG,
                      dplyr::matches('pct_tg\\d+|^tg\\d+')) %>%
        dplyr::mutate_each(dplyr::funs(as.numeric(scale(.))), -SID)

    gee_ready_data <-
        dplyr::full_join(
            no_fattyacids,
            scaled_variables,
            by = 'SID'
        ) %>%
        dplyr::group_by(VN) %>%
        dplyr::mutate(
            Waist = as.numeric(scale(Waist)),
            ALT = as.numeric(scale(ALT)),
            BaseAge = as.numeric(scale(BaseAge)),
            MET = as.numeric(scale(MET))
        ) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(SID, VN)

    return(gee_ready_data)
}

prep_dysglycemia_data <- function(data) {
    dysgly.data <-
        dplyr::left_join(
            data %>%
                dplyr::filter(VN == 0),
            data %>%
                dplyr::filter(!is.na(TotalNE)) %>%
                dplyr::mutate_each(dplyr::funs(ifelse(is.na(.), 0, .)), IFG, IGT) %>%
                dplyr::mutate(PreDM = as.numeric(rowSums(.[c('IFG', 'IGT')], na.rm = TRUE))) %>%
                dplyr::mutate(FactorDysgly = ifelse(
                    PreDM == 1, 'PreDM',
                    ifelse(DM == 1, 'DM',
                           'NGT')
                )) %>%
                dplyr::select(SID, VN, FactorDysgly) %>%
                tidyr::spread(VN, FactorDysgly) %>%
                dplyr::mutate(
                    DetailedConvert = as.factor(paste0(`0`, '-', `1`, '-', `2`)),
                    ConvertDysgly = as.numeric(!grepl('NGT$|NGT-NA$|NGT-NGT-NA$|NGT-NA-NA$',
                                                      DetailedConvert)),
                    ConvertDM = as.numeric(grepl('-DM$|-DM-DM$|-DM-NA$|-DM-NGT$',
                                                 DetailedConvert)),
                    ConvertPreDM = as.numeric(grepl('-PreDM$|-PreDM-PreDM$|-PreDM-NA$',
                                                    DetailedConvert))
                )
        ) %>%
        dplyr::filter(!is.na(TotalNE))

    return(dysgly.data)
}
