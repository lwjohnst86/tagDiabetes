tg_conc <- grep('^tg\\d+', vars, value = TRUE)
tg_pct <- grep('^pct_tg\\d+', vars, value = TRUE)
tg_totals <- c('TotalTG', 'BaseTAG')
outcomes <- c('lHOMA2_S', 'lISI', 'lIGIIR', 'lISSI2')
covariates <- c('YearsFromBaseline', 'Waist', 'BaseAge', 'BiEthnicity', 'Sex', 'ALT', 'MET', 'TotalNE')
