# Plot --------------------------------------------------------------------

#' Interaction by time plots using GEE.
#'
#' @param data Project data
#' @param y Outcome variable.
#' @param x Predictor variable (Fatty acid).
#' @param fatty.acid.levels Whether to split into median, tertiles, or quartiles.
#' @param covars Covariates to adjust for.
#' @param title The title letter of the plot.
#' @param interaction.term Which variable to use for the interaction (VN).
#' @param is.fatty.acid Whether the x is a fatty acids.
#' @param y.axis.label The y axis label.
#'
#' @export
plot_interaction <- function(data, y, x,
                         fatty.acid.levels = c(2, 3, 4),
                         covars = covariates,
                         title = 'A',
                         interaction.term = 'VN',
                         is.fatty.acid = TRUE,
                         y.axis.label = NULL) {
    fa.lvl <- fatty.acid.levels
    if (!fa.lvl %in% c(2, 3, 4)) {
        print('Please select either 2, 3, or 4 for the fatty.acid.levels')
    } else if (fa.lvl == 2) {
        lvls <- 0:2/2
        lvl.labs <- c('<50%', '>50%')
    } else if (fa.lvl == 3) {
        lvls <- 0:3/3
        lvl.labs <- c('<33%', '33-66%', '>66%')
    } else if (fa.lvl == 4) {
        lvls <- 0:4/4
        lvl.labs <- c('Lowest', 'Low', 'High', 'Highest')
    }

    if (is.null(y.axis.label)) {
        y.axis.label <- paste0('log(', substring(y, 2), ')')
    }

    fatty.acid.name <- renaming_fats(x) %>%
      gsub('(\\D\\D)\\.', '\\U\\1 ', ., perl = TRUE)

    if (is.fatty.acid) {
        fatty.acid.name <- paste(fatty.acid.name,
                                 ifelse(grepl('pct_', x), '(mol%)', '\n(nmol/mL)'))
    }

    geeFormula <- paste(y, 'f.fat', sep = ' ~ ') %>%
      paste(., paste(covars, collapse = ' + '),
            paste('f.fat', interaction.term, sep = ':'), sep = ' + ') %>%
      as.formula()

    prep.fit <- data %>%
      dplyr::select_(.dots = c(y, x, covars, interaction.term, 'SID')) %>%
      dplyr::filter(complete.cases(.)) %>%
      dplyr::rename_('fa' = x, 'VN' = interaction.term) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(f.fat = cut(fa, stats::quantile(fa, probs = lvls),
                         include.lowest=TRUE, labels = lvl.labs),
             VN = factor(VN, ordered = TRUE))

    fit <-
        geepack::geeglm(geeFormula, id = SID, family = stats::gaussian,
               corstr = 'ar1', data = prep.fit)

    dodge <- ggplot2::position_dodge(width = 0.1)

    lsmeans::lsmeans(fit, ~ f.fat:VN) %>%
      summary() %>%
      as.data.frame() %>%
      ggplot2::ggplot(ggplot2::aes(VN, lsmean, group = f.fat, shape = f.fat, ymax = max(asymp.UCL))) +
      ggplot2::geom_line(ggplot2::aes(linetype = f.fat), position = dodge) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = asymp.LCL, ymax = asymp.UCL, linetype = f.fat), width = 0.1,
                    position = dodge) +
      ggplot2::geom_point(position = dodge) +
      ggplot2::scale_shape_discrete('Percentile') +
      ggplot2::scale_linetype_discrete('Percentile') +
      ggplot2::scale_x_discrete('Visit number', labels = c('0-yr', '3-yrs', '6-yrs')) +
      ggthemes::theme_tufte(10, base_family = 'Arial') +
      ggplot2::theme(legend.position = 'bottom',
            strip.background = ggplot2::element_rect(fill = 'grey95', colour = 'grey95'),
            plot.margin = grid::unit(c(0.5, 0, 0, 0), "cm"),
            legend.key.width = grid::unit(0.75, "line"),
            legend.key.height = grid::unit(0.75, "line")
            ) +
      ## theme(legend.justification = c(1, 1),
      ##       legend.position = c(1, 1)) +
      ggplot2::ylab(y.axis.label) +
      ggplot2::ggtitle(paste0(title, ': ', fatty.acid.name)) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0, vjust = 1))
}
