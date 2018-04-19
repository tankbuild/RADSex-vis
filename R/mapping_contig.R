#' @title Mapping contig
#'
#' @description Generates a linear plot of RADSex mapping results for a specified contig. The resulting figure contains two plot: the top plot shows sex-bias against
#' position on the contig, and the bottom plot shows probability of association with sex against position on the contig.
#' A specific region of the contig can also be plotted.
#'
#' @param data A list of mapping results obtained with the \code{\link{load_mapping_results}} function.
#'
#' @param contig Name of the contig to plot.
#'
#' @param region A vector of two integers specifying the region of the contig to plot (default NULL).
#'
#' @param signif.threshold Significance threshold for association with sex (default 0.05).
#'
#' @param point.size Size of a point in the plot (default 1.5).
#'
#' @param color.sex.bias If TRUE, points on the sex-bias track will be colored according to sex.bias.palette (default TRUE).
#'
#' @param sex.bias.palette A vector of three colors defining the sex-bias track palette: female-biased, neutral, male-biased. (default c("firebrick1", "black", "dodgerblue2"))
#'
#' @examples
#'
#' contig_lengths <- load_contig_lengths("contig_lengths.tsv")
#' contig_names <- load_contig_names("contig_names.tsv")
#' data <- load_mapping_results("mapping_results.tsv", contig_lengths, contig_names = contig_names, plot.unplaced = FALSE)
#'
#' contig_plot <- mapping_contig(data, "LG05", region = c(1500000, 2500000), point.size = 1, color.sex.bias = FALSE)


mapping_contig <- function(data, contig, region = NULL,
                           signif.threshold = 0.05,
                           point.size = 1.5,
                           color.sex.bias = TRUE, sex.bias.palette = c("firebrick1", "grey10", "dodgerblue2")) {


    # Check if the contig exists in the mapping results
    if (contig %in% names(data$lengths)) {

        contig_name <- contig  # If the contig name was specified as the ID in the reference genome, save it as it is

    } else if (!is.null(data$names) & contig %in% data$names) {

        contig_name <- names(data$names)[which(data$names == contig)]  # If the contig name was specified as one of the chromosomes names, save the corresponding ID

    } else {

        stop(paste0("Could not find contig \"", contig , "\"."))

    }

    # Adjust significance threshold for bonferroni correction
    signif.threshold <- -log(signif.threshold / dim(data$data)[1], 10)

    # Get contig info
    x_limits <- c(0, data$lengths[contig_name])
    contig_data <- subset(data$data, data$data$Contig == contig_name)

    # Adjust contig info if a region was specified
    if (!is.null(region)) {
        x_limits <- region
        contig_data <- subset(contig_data, contig_data$Position >= region[1] & contig_data$Position <= region[2])
    }

    # Generate a good x axis scale according to the size of the region
    x_scale <- generate_x_scale(x_limits, contig)

    # Generate the top plot (sex bias)
    g <- ggplot2::ggplot(contig_data, ggplot2::aes(x = Position, y = SexBias, color = SexBias)) +
        cowplot::theme_cowplot() +
        ggplot2::geom_point(size = point.size) +
        ggplot2::scale_y_continuous(name = "Sex Bias", limits = c(-1, 1)) +
        ggplot2::theme(axis.title.x = ggplot2::element_blank(), legend.position = "none") +
        x_scale

    # Add color scale if specified
    if (color.sex.bias) {
        g <- g + ggplot2::scale_color_gradientn(name = "Sex Bias", colours = colorRampPalette(sex.bias.palette)(20),
                                                limits = c(-1, 1))
    }

    # Generate the bottome plot (probability of association with sex)
    h <- ggplot2::ggplot(contig_data, ggplot2::aes(x = Position, y = P)) +
        cowplot::theme_cowplot() +
        ggplot2::geom_point(size = point.size, color = "grey20") +
        ggplot2::scale_y_continuous(name = expression(paste0('-log'['10'], '(p'['Association with sex'], ')')),
                                    limits = c(0, max(contig_data$P, 1.25 * signif.threshold))) +
        ggplot2::geom_hline(yintercept = signif.threshold, color = "red3", lty = 1, lwd = 1) +
        ggplot2::annotate("text", x = x_limits[1] + 0.05 * (x_limits[2] - x_limits[1]), y = signif.threshold + 0.075 * max(contig_data$P, 1.25 * signif.threshold),
                          label = "p = 0.05", color = "red3", size = 5) +
        x_scale

    # Combine the two plots
    combined <- cowplot::plot_grid(g, h, ncol = 1, align = "v")

    return(combined)
}


# Generate a nice scale for any interval
generate_x_scale <- function(region, contig_name) {

    # Size of the region
    S <- (region[2] - region[1]) / 10
    if (S <= 0) stop(paste0("Error: the size of the region has to be > 0 (value: ", S * 10, ")"))

    # Find the order of magnitude of the region's size
    N <- floor(log(S, 10))
    N10 <- 10 ^ N

    # Generate a scale of 10 round values within the region
    scale <- seq(N10 * floor(region[1] / N10), N10 * ceiling(region[2] / N10), N10 * round(S / N10, 1))

    # Adjust the labels based on the size of the values (megabp or kilop)
    if (region[1] < 10 ^ 6 & region[2] < 10 ^ 6) {
        scale_labels <- round(scale / 10 ^ 3, 3 - N)
        bp_unit <- "K"
    } else {
        scale_labels <- round(scale / 10 ^ 6, 6 - N)
        bp_unit <- "M"
    }

    output <- ggplot2::scale_x_continuous(name = paste0("Position on ", contig_name, " (", bp_unit, "bp.)"),
                                          breaks = scale, labels = scale_labels, limits = c(min(scale[1], region[1]), max(tail(scale, 1), region[2])))

    return(output)
}

