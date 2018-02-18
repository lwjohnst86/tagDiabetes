# Renaming ----------------------------------------------------------------

renaming_table_rows <- function(x) {
    x %>%
        gsub('IGIIR', 'IGI/IR', .) %>%
        gsub('ISSI2', 'ISSI-2', .) %>%
        gsub('HOMA2IR', 'HOMA2-IR', .) %>%
        gsub('HOMA2_S', 'HOMA2-%S', .) %>%
        gsub('TAG', 'TG (mmol/L)', .) %>%
        gsub('Chol', 'Chol (mmol/L)', .) %>%
        gsub('LDL', 'LDL (mmol/L)', .) %>%
        gsub('ALT', 'ALT (U/L)', .) %>%
        gsub('HDL', 'HDL (mmol/L)', .) %>%
        gsub('BaseTotalTG', 'TGFA (nmol/mL)', .) %>%
        gsub('BaseTotalNE', 'NEFA (nmol/mL)', .) %>%
        gsub('Age', 'Age (yrs)', .) %>%
        gsub('BMI', 'BMI (kg/m^2^)', .) %>%
        gsub('Waist', 'WC (cm)', .)
}

renaming_outcomes <- function(x) {
    x %>%
        gsub('linvHOMA', 'log(1/HOMA-IR)', .) %>%
        gsub('lHOMA2IR', 'log(HOMA2-IR)', .) %>%
        gsub('lHOMA2_S', 'log(HOMA2-%S)', .) %>%
        gsub('lISI', 'log(ISI)', .) %>%
        gsub('lIGIIR', 'log(IGI/IR)', .) %>%
        gsub('lISSI2', 'log(ISSI-2)', .) %>%
        gsub('^HOMA2_S$', 'HOMA2-%S', .) %>%
        gsub('^ISSI2$', 'ISSI-2', .) %>%
        gsub('^IGIIR$', 'IGI/IR', .)
}

renaming_fats <- function(x) {
    x %>%
        stringr::str_replace('^pct_', '') %>%
        stringr::str_replace('TotalTG', 'Total') %>%
        stringr::str_replace('BaseTAG', 'Clinical TG') %>%
        PROMISE.misc::renaming_fa(keep.fraction = FALSE)
}

renaming_fraction <- function(x) {
    x %>%
        gsub('tg', 'Triacylglycerol', .)
}
