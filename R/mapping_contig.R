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
#' @param title Title of the plot (default NULL).
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
#' data <- load_mapping_results("mapping_results.tsv", contig_lengths_file_path = 'contig_lengths.tsv',
#'                              contig_names_file_path = 'contig_names.tsv')
#'#' contig_plot <- mapping_contig(data, "LG05", region = c(1500000, 2500000), point.size = 1, color.sex.bias = FALSE)


mapping_contig <- function(data, contig, region = NULL, title = NULL,
                           signif.threshold = 0.05, point.size = 1.5,
                           color.sex.bias = TRUE, sex.bias.palette = c("firebrick1", "grey10", "dodgerblue2")) {

    # Getting contig information

    if (contig %in% data$names) {

        # Case of contig = "LG05"
        contig_name <- contig
        contig <- names(data$names[which(data$names == contig)])

    } else if (contig %in% names(data$names)) {

        # Case of contig = "NC_0215354.2"
        contig_name <- data$names[which(names(data$names) == contig)]

    } else if (contig %in% c(names(data$lengths$unplaced), names(data$lengths$lg))) {

        # Case of unplaced scaffold or no contig names
        contig_name <- contig

    } else {

        stop(paste0(" - Error: contig \"", contig, "\" does not exist."))

    }

    # Adjust significance threshold for bonferroni correction
    signif.threshold <- -log(signif.threshold / dim(data$data)[1], 10)

    # Setting region to entire contig if region is NULL
    if (is.null(region)) {
        if (contig %in% names(data$lengths$lg)) {
            region <- c(0, data$lengths$lg[[contig]])
        } else if (contig %in% names(data$lengths$unplaced)){
            region <- c(0, data$lengths$unplaced[[contig]])
        } else {
            stop(paste0(" - Error: could not find length for contig \"", contig, "\"."))
        }
    }

    # Extract region from contig data
    contig_data <- subset(data$data, data$data$Contig_id == contig &
                              data$data$Original_position >= region[1] &
                              data$data$Original_position <= region[2])

    # Generate a good x axis scale according to the size of the region
    x_scale <- generate_x_scale(region, contig_name)

    # Generate the top plot (sex bias)
    g <- ggplot2::ggplot(contig_data, ggplot2::aes(x = Original_position, y = SexBias, color = SexBias)) +
        cowplot::theme_cowplot() +
        ggplot2::geom_point(size = point.size) +
        ggplot2::scale_y_continuous(name = "Sex Bias", limits = c(-1, 1)) +
        ggplot2::theme(axis.title.x = ggplot2::element_blank(), legend.position = "none") +
        ggplot2::geom_hline(yintercept = 0, linetype = 3, color = "grey50") +
        x_scale

    # Add color scale if specified
    if (color.sex.bias) {
        g <- g + ggplot2::scale_color_gradientn(name = "Sex Bias", colours = colorRampPalette(sex.bias.palette)(20),
                                                limits = c(-1, 1))
    }

    # Add title if specified
    if (!is.null(title)) {
        g <- g + ggplot2::ggtitle(title) +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 20, face = "bold", margin = ggplot2::margin(0, 0, 10, 0)))
    } else {
        g <- g + ggplot2::theme(plot.title = ggplot2::element_blank())
    }

    # Generate the bottome plot (probability of association with sex)
    h <- ggplot2::ggplot(contig_data, ggplot2::aes(x = Original_position, y = P)) +
        cowplot::theme_cowplot() +
        ggplot2::geom_point(size = point.size, color = "grey20") +
        ggplot2::scale_y_continuous(name = expression(paste('-log'['10'], '(p'['Association with sex'], ')', sep = "")),
                                    limits = c(0, max(contig_data$P, 1.25 * signif.threshold))) +
        ggplot2::geom_hline(yintercept = signif.threshold, color = "red3", lty = 1, lwd = 1) +
        ggplot2::annotate("text", x = region[1] + 0.05 * (region[2] - region[1]), y = signif.threshold + 0.075 * max(contig_data$P, 1.25 * signif.threshold),
                          label = "p = 0.05", color = "red3", size = 5) +
        x_scale

    # Combine the two plots
    if (!is.null(title)) {

        combined <- cowplot::plot_grid(g, h, ncol = 1, align = "v", rel_heights = c(1.1, 1))

    } else {

        combined <- cowplot::plot_grid(g, h, ncol = 1, align = "v")

    }

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

    output <- ggplot2::scale_x_continuous(name = paste0("Position on ", contig_name, " (", bp_unit, "bp)"),
                                          breaks = scale, labels = scale_labels, limits = c(min(scale[1], region[1]), max(tail(scale, 1), region[2])))

    return(output)
}

