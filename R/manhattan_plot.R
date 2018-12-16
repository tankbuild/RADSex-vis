#' @title Mapping manhattan plot
#'
#' @description Generates a manhattan plot of association with sex from RADSex mapping results.
#'
#' @param data A list of mapping results obtained with the \code{\link{load_mapping_results}} function.
#'
#' @param point.size Size of a point in the plot (default 0.5)
#'
#' @param point.palette Color palette for the dots (default c("dodgerblue3", "darkgoldenrod2"))
#'
#' @param background.palette Color palette for the background (default c("grey85", "grey100"))
#'
#' @param signif.threshold Significance threshold for association with sex (default 0.05).
#'
#' @param significance.line.color Color for significance line, set to NULL for no line (default "black").
#'
#' @param significance.line.type Linetype for the significance line, as usually defined in R (default 2).
#'
#' @param significance.text.position X and Y axis offset for the significance text, as fractions of total axis length (default c(0.05, 0.05)).
#'
#' @examples
#'
#' data <- load_mapping_results("mapping_results.tsv", "contig_lengths.tsv", contig_names = "contig_names.tsv", plot.unplaced = FALSE)
#'
#'  mapping_manhattan_plot(data, point.size = 0.5, signif.threshold = 0.05,
#'                         point.palette = c("dodgerblue3", "darkgoldenrod2"),
#'                         background.palette = c("grey85", "grey100"),
#'                         significance.line.color = "black",
#'                         significance.line.type = 2,
#'                         significance.text.position = c(0.05, 0.05)


mapping_manhattan_plot <- function(data, point.size = 0.5, signif.threshold = 0.05,
                                   point.palette = c("dodgerblue3", "darkgoldenrod2"),
                                   background.palette = c("grey85", "grey100"),
                                   significance.line.color = "black",
                                   significance.line.type = 2,
                                   significance.text.position = c(0.05, 0.05)) {

    # Apply bonferroni correction for significance threshold
    signif_threshold <- -log(signif.threshold / dim(data$data)[1], 10)

    # Compute cumulative lengths of contigs, which will be added to the position of each point depending on the contig
    cumulative_lengths <- c(0, cumsum(data$lengths$lg))
    names(cumulative_lengths) <- c(names(data$lengths$lg), "Unplaced")

    # Adjust x-axis position for each point in the data based on cumulative lengths of contigs
    manhattan_data <- data$data
    manhattan_data$Position = manhattan_data$Position + cumulative_lengths[manhattan_data$Contig]

    # Attribute alternating colors to each contig
    order <- seq(1, length(data$lengths$plot))
    names(order) <- names(data$lengths$plot)
    manhattan_data$Color <- order[manhattan_data$Contig] %% 2
    manhattan_data$Color <- as.factor(as.character(manhattan_data$Color))

    # Create the background data for alternating background color
    background = data.frame(start = cumulative_lengths, end = cumulative_lengths + data$lengths$plot)
    background$Color = rep_len(c("A", "B"), length.out = dim(background)[1])  # Alternating A/B

    # Create merged color palette for background and data points
    merged_color_palette <- c("0"=point.palette[1], "1"=point.palette[2], "B"=background.palette[1], "A"=background.palette[2])

    ymax = 1.1 * max(manhattan_data$P, 1.1 * signif_threshold)

    manhattan_plot <- ggplot2::ggplot() +
        # Backgrounds with alternating colors
        ggplot2::geom_rect(data = background,
                           ggplot2::aes(xmin = start, xmax = end, ymin = 0, ymax = ymax, fill = Color),
                           alpha = 0.5) +
        # Data points
        ggplot2::geom_point(data = manhattan_data,
                            ggplot2::aes(x = Position, y = P, color = Color),
                            size = point.size,
                            alpha = 1) +
        # Attribute color values from merged color scale for points and backgrounds
        ggplot2::scale_color_manual(values = merged_color_palette) +
        ggplot2::scale_fill_manual(values = merged_color_palette) +
        # Generate x-axis, use background start and end to place LG labels
        ggplot2::scale_x_continuous(name="Linkage Group",
                                    breaks = background$start + (background$end - background$start) / 2,
                                    labels = names(order),
                                    expand = c(0, 0)) +
        # Generate y-axis
        ggplot2::scale_y_continuous(name=expression(paste('-log'['10'], '(p'['Association with sex'], ')', sep = "")),
                                    limits=c(-0.05, ymax),
                                    expand = c(0, 0)) +
        # Adjust theme elements
        cowplot::theme_cowplot() +
        ggplot2::theme(legend.position = "none",
                       panel.grid.minor = ggplot2::element_blank(),
                       panel.grid.major.x = ggplot2::element_blank(),
                       axis.line.x = ggplot2::element_blank(),
                       panel.border = ggplot2::element_blank(),
                       axis.line.y = ggplot2::element_line(color = "black"),
                       axis.title.x = ggplot2::element_text(face="bold", margin = ggplot2::margin(10, 0, 0, 0)),
                       axis.title.y = ggplot2::element_text(face="bold", margin = ggplot2::margin(0, 10, 0, 0)),
                       axis.text = ggplot2::element_text(color="black", face="bold")) +
        # Draw line for significance threshold
        ggplot2::geom_hline(yintercept = signif_threshold, color = significance.line.color, lty = significance.line.type, lwd = 0.5) +
        ggplot2::annotate("text", x = significance.text.position[1] * max(manhattan_data$Position),
                          y = signif_threshold + significance.text.position[2] * ymax,
                          label = paste("p = ", as.character(signif.threshold)), color = significance.line.color, size = 5)

    return(manhattan_plot)
}
