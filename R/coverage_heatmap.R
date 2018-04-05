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
#' @param individual.names If TRUE, shows individual names on the x-axis.
#'
#' @param sequence.names If TRUE, shows sequence names on the y-axis.
#'
#' @param individual.dendrogram If TRUE, shows individual clustering dendrogram on the x-axis.
#'
#' @param sequence.dendrogram If TRUE, shows sequence clustering dendrogram on the y-axis.
#'
#' @return A coverage heatmap stored in a ggplot object.
#'
#' @examples
#' heatmap = coverage_heatmap(data, popmap=popmap,
#'                            title="Individuals and sequences clustering based on coverage",
#'                            males.color="blue", females.color="red",
#'                            coverage.palette=c("white", "green", "black", "red", "red3"),
#'                            individual.names=TRUE, sequence.names=TRUE,
#'                            individual.dendrogram=TRUE, sequence.dendrogram=TRUE)


coverage_heatmap <- function(data, popmap=NULL, title=NULL,
                             males.color="dodgerblue3", females.color="red3",
                             coverage.palette=c("white", "royalblue2", "black", "gold2", "red3"),
                             individual.names=TRUE, sequence.names=FALSE,
                             individual.dendrogram=TRUE, sequence.dendrogram=TRUE) {

    # Define the color palette for individual names: black if not popmap, or colored by sex if popmap
    individual_names = data$individuals$labels[data$individuals$order]
    if (is.null(popmap)) {
        sex_palette = rep("black", length(individual_names))
        names(sex_palette) = individual_names
    } else {
        temp = c("M"=males.color, "F"=females.color)
        sex_palette = temp[popmap[individual_names]]
    }

    # Compute the main heatmap object
    heatmap <- ggplot2::ggplot(data$data, ggplot2::aes(x = individual, y = id, fill = coverage)) +
        ggplot2::geom_tile(color = "grey30", size = 0.02) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_blank(),
              axis.text.y = ggplot2::element_blank(),
              axis.ticks.y = ggplot2::element_blank(),
              axis.title.y = ggplot2::element_blank(),
              axis.title.x = ggplot2::element_blank(),
              plot.margin = ggplot2::margin(15, 15, 15, 30),
              panel.border = ggplot2::element_rect(size = 0.75, color = "black"),
              legend.position = "right",
              legend.background = ggplot2::element_rect(color="black", size=0.5),
              legend.margin = ggplot2::margin(5, 5, 5, 5),
              legend.title = ggplot2::element_text(size = 14, face = "bold"),
              legend.text = ggplot2::element_text(size = 11),
              legend.key.height = ggplot2::unit(0.06, "npc"),
              legend.key.width = ggplot2::unit(0.04, "npc")) +
        ggplot2::scale_fill_gradientn(name = "Cov.", colours = coverage.palette, values = c(0, 0.00001, data$distribution[4]/data$distribution[6], data$distribution[5]/data$distribution[6], 1)) +
        ggplot2::scale_x_discrete(expand = c(0, 0)) +
        ggplot2::scale_y_discrete(expand = c(0, 0))

    # Add title if specified
    if (!is.null(title)) {
        heatmap <- heatmap + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 20, face = "bold", margin = ggplot2::margin(0, 0, 10, 0))) +
            ggplot2::ggtitle(title)
    } else {
        heatmap <- heatmap + ggplot2::theme(plot.title = ggplot2::element_blank())
    }

    # Base position of dendrograms in the gtable
    individual_dendrogram_top <- 7
    sequence_dendrogram_top <- 4

    # Add individual names if specified
    if (individual.names) {
        heatmap <- heatmap + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, colour = sex_palette, size = 10, margin = ggplot2::margin(2.5,0,0,0)))
        individual_dendrogram_top <- individual_dendrogram_top + 1  # Increment position of individual dendrogram in the gtable because one row was added.
        }

    # Add sequence names if specified
    if (sequence.names) {
        heatmap <- heatmap + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10))
    }

    # Create plot grob of heatmap for the combined gtable
    combined <- ggplot2::ggplotGrob(heatmap)

    if (individual.dendrogram == TRUE) {

        # Compute individual dendrogram object
        individual_dendrogram <- suppressMessages(ggdendro::ggdendrogram(individual_clusters, labels = FALSE, leaf_labels = FALSE, theme_dendro = TRUE, rotate = FALSE) +
                                                      ggplot2::theme(plot.margin = grid::unit(c(0.1, 0.01, 0, 0.01), 'npc'),
                                                                     axis.text.x = ggplot2::element_blank(),
                                                                     axis.text.y = ggplot2::element_blank(),
                                                                     axis.title.x = ggplot2::element_blank()) +
                                                      ggplot2::scale_y_reverse(expand = c(0, 0.5)) +
                                                      ggplot2::scale_x_continuous(expand = c(0, 0)))

        # Add a row to the combined gtable for the dendrogram
        combined <- gtable::gtable_add_rows(combined, grid::unit(0.1, "npc"), pos=individual_dendrogram_top)
        # Add the dendrogram to the combined gtable
        combined <- gtable::gtable_add_grob(combined, ggplot2::ggplotGrob(individual_dendrogram), t = individual_dendrogram_top, l=4, b=individual_dendrogram_top + 1, r=5)

        sequence_dendrogram_top <- sequence_dendrogram_top + 1
    }

    if (sequence.dendrogram) {

        sequence_dendrogram <- suppressMessages(ggdendro::ggdendrogram(sequence_clusters, labels = FALSE, leaf_labels = FALSE, theme_dendro = TRUE, rotate = FALSE) +
                                                    ggplot2::theme(plot.margin = grid::unit(c(0.05, 0.1, 0.05, 0), 'npc'),
                                                                   axis.text.x = ggplot2::element_blank(),
                                                                   axis.text.y = ggplot2::element_blank(),
                                                                   axis.title.x = ggplot2::element_blank()) +
                                                    ggplot2::coord_flip() +
                                                    ggplot2::scale_y_reverse(expand = c(0, 0.5)) +
                                                    ggplot2::scale_x_continuous(expand=c(0, 0)))

        # Add a row to the combined gtable for the dendrogram
        combined <- gtable::gtable_add_cols(combined, grid::unit(0.04, "npc"), pos=0)
        # Add the dendrogram to the combined gtable
        combined <- gtable::gtable_add_grob(combined, ggplot2::ggplotGrob(sequence_dendrogram), t = sequence_dendrogram_top, l=1, b=sequence_dendrogram_top + 1, r=2)
    }

    return(combined)
}
