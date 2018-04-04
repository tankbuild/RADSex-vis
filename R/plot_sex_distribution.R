#' @title Plot sex distribution
#'
#' @description Plot the heatmap of distribution of sequences between sexes from a table generated with RADSex.
#'
#' @param input_file_path Path to a table of distribution of sequences between sexes.
#'
#' @param output_file_path Path to the plot output file. If no output file is specified, the heatmap will be plotted in the default device.
#'
#' @param width Width of the output file (in inches).
#'
#' @param height Height of the output file (in inches).
#'
#' @param dpi Resolution of the output file.
#'
#' @param autoscale If set to TRUE, the height of the output file will be automatically set so that horizontal and vertical scales are equal.
#' The resulting heatmap will not be square if the number of males and females in the population are very different.
#'
#' @param title Plot title.
#'
#' @param show.significance If TRUE, tiles with significant association with sex are highlighted with the color defined in the significance.color parameter.
#'
#' @param significance.color Color of the border for tiles significantly associated with sex.
#'
#' @param significance.threshold P-value threshold to consider a tile significantly associated with sex.
#'
#' @param color.scale.bins A vector of values to use as bins in the color palette.
#'
#' @param color.scale.colors A vector of two colors used to create the color palette gradient.
#'
#' @return A heatmap stored in a ggplot object.
#'
#' @examples
#' # No autoscale: the plot dimensions are specified with "width" and "height".
#' heatmap = sex_distribution_heatmap(input_file_path, output_file_path='heatmap.png',
#'                                    title="Distribution of sequences between sexes",
#'                                    width=10, height=10, dpi=200,
#'                                    significance.color="gold", significance.threshold=0.05,
#'                                    color.scale.bins=c(0, 10, 100, 1000),
#'                                    color.scale.colors=c("white", "grey10"))
#'
#' # Autoscale: the plot dimensions are automatically calculated from "width".
#' heatmap = sex_distribution_heatmap(input_file_path, output_file_path='heatmap.png',
#'                                    title="Distribution of sequences between sexes",
#'                                    width=10, autoscale=TRUE)

plot_sex_distribution <- function(input_file_path, output_file_path=NULL, title=NULL,
                                  width=8, height=6, dpi=300, autoscale=FALSE,
                                  show.significance=TRUE, significance.color="red3", significance.threshold=0.05,
                                  color.scale.bins=c(0, 1, 5, 25, 100, 1000), color.scale.colors=c("white", "navyblue")) {

    if (!file.exists(input_file_path)) {
        stop(paste0("The specified input file (", input_file_path, ") does not exist."))
    }

    data <- load_sex_distribution_table(input_file_path)
    heatmap <- sex_distribution_heatmap(data,  title=title,
                                        show.significance=show.significance, significance.color=significance.color, significance.threshold=significance.threshold,
                                        color.scale.bins=color.scale.bins, color.scale.colors=color.scale.colors)

    if (autoscale) {
        ratio <- data$n_females / data$n_males
        if (is.null(title)) {
            height <- width * (ratio - 0.05)
        } else {
            height <- width * ratio
        }
    }

    if (!is.null(output_file_path)) {
        ggplot2::ggsave(output_file_path, plot=heatmap, width=width, height=height, dpi=dpi)
    } else {
        print(heatmap)
    }
}
