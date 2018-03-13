# Misc --------------------------------------------------------------------

graph_theme <- function(base.plot, ticks = TRUE, minor.grid.lines = FALSE, legend.pos = 'bottom') {
    graph.theme <-
        "%+replace%"(
            ggthemes::theme_tufte(base_size = 10, base_family = 'sans'),
            theme(
                axis.line = element_line('black'),
                axis.line.x = element_line('black'),
                axis.line.y = element_line('black'),
                legend.key.width = grid::unit(0.7, "line"),
                legend.key.height = grid::unit(0.7, "line"),
                strip.background = element_blank(),
                plot.margin = grid::unit(c(0.5, 0, 0, 0), "cm"),
                legend.position = legend.pos,
                strip.placement = "outside"
            )
        )

    if (!ticks) {
        graph.theme <- "%+replace%"(graph.theme,
                                             theme(axis.ticks.y = element_blank()))
    }

    if (minor.grid.lines) {
        graph.theme <- "%+replace%"(
            graph.theme,
            theme(
                panel.grid = element_line(),
                panel.grid.minor = element_blank()
            )
        )
    }

    return(graph.theme)
}
