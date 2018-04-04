#' @title Coverage heatmap
#'
#' @description Generates a heatmap of coverage for a subset of sequences. In the resulting heatmap, the color of a tile at coordinates (x, y)
#' indicates the coverage of the sequence y in the individual x.
#'
#' @param data A table of coverage obtained with the \code{\link{load_coverage_table}} function.
#'
#' @param popmap A population map obtained with the \code{\link{load_population_map}} function.
#'
#' @param title Plot title.
#'
#' @param males.color If a popmap is specified, sets the color of male individual names on the horizontal axis.
#'
#' @param females.color If a popmap is specified, sets the color of female individual names on the horizontal axis.
#'
#' @param coverage.palette Color palette for coverage. The value should be a vector of length 5 corresponding to the following intervals/values:
#' (0, 1 : mean, mean : 3rd quartile, 3rd quartile : (max - 1), max)
#'
#' @return A coverage heatmap stored in a ggplot object.
#'
#' @examples
#' heatmap = coverage_heatmap(data, popmap=popmap,
#'                            title="Individuals and sequences clustering based on coverage",
#'                            males.color="blue", females.color="red",
#'                            coverage.palette=c("white", "green", "black", "red", "red3"))

coverage_heatmap <- function(data, popmap=NULL, title=NULL,
                             males.color="dodgerblue3", females.color="red3",
                             coverage.palette=c("white", "royalblue2", "black", "gold2", "red3")) {

    individual_names = data$individuals$labels[data$individuals$order]
    if (is.null(popmap)) {
        sex_palette = rep("black", length(individual_names))
        names(sex_palette) = individual_names
    } else {
        temp = c("M"=males.color, "F"=females.color)
        sex_palette = temp[popmap[individual_names]]
    }

    heatmap <- ggplot2::ggplot(data$data, ggplot2::aes(x = individual, y = id, fill = coverage)) +
        ggplot2::geom_tile(color = "grey30", size = 0.02) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, colour = sex_palette, size = 10),
              axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(),
              axis.title.x = ggplot2::element_blank(), plot.margin = ggplot2::margin(15, 15, 15, 30),
              panel.border = ggplot2::element_rect(size = 0.75, color = "black"),
              legend.position = "right", legend.margin = ggplot2::margin(0, 0, 0, 0),
              legend.title = ggplot2::element_text(size = 14, face = "bold"), legend.text = ggplot2::element_text(size = 11),
              legend.key.height = ggplot2::unit(0.06, "npc"), legend.key.width = ggplot2::unit(0.04, "npc")) +
        ggplot2::scale_fill_gradientn(name = "Cov.", colours = coverage.palette, values = c(0, 0.00001, data$distribution[4]/data$distribution[6], data$distribution[5]/data$distribution[6], 1))

    if (!is.null(title)) {
        heatmap <- heatmap + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 20, face = "bold", margin = ggplot2::margin(0, 0, 10, 0))) +
            ggplot2::ggtitle(title)
    } else {
        heatmap <- heatmap + ggplot2::theme(plot.title = ggplot2::element_blank())
    }

    combined <- ggplot2::ggplotGrob(heatmap)

    individual_dendrogram <- suppressMessages(ggdendro::ggdendrogram(individual_clusters, labels = FALSE, leaf_labels = FALSE, theme_dendro = TRUE, rotate = FALSE) +
                                                  ggplot2::theme(plot.margin = ggplot2::unit(c(0.1, 0, 0, 0), 'cm'), axis.text.x = ggplot2::element_blank(),
                                                                 axis.text.y = ggplot2::element_blank(), axis.title.x = ggplot2::element_blank()) +
                                                  ggplot2::scale_y_reverse(expand = c(0,0.5)) + ggplot2::scale_x_continuous(expand = c(0,0)))

    combined <- gtable::gtable_add_rows(combined, grid::unit(0.06, "npc"), pos=8)
    combined <- gtable::gtable_add_grob(combined, ggplot2::ggplotGrob(individual_dendrogram), t = 8, l=4, b=9, r=5)

    return(combined)
}
