# Grab or combine data ----------------------------------------------------

#' Wrangle data to determine dysglycemia conversion.
#'
#' @param data The project data.
#'
#' @export
prep_dys_data <- function(data) {
    dysgly.data <-
        dplyr::left_join(
            data %>%
                dplyr::filter(VN == 0),
            data %>%
                dplyr::filter(!is.na(TotalTG)) %>%
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
        dplyr::filter(!is.na(TotalTG))

    return(dysgly.data)
}

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
