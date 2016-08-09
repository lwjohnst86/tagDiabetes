tg_conc <- grep('^tg\\d+', names(project_data), value = TRUE)
tg_pct <- grep('^pct_tg\\d+', names(project_data), value = TRUE)
tg_totals <- c('TotalTG', 'TAG')
outcomes <- c('linvHOMA', 'lISI', 'lIGIIR', 'lISSI2')
