#' @title Sex distribution heatmap
#'
#' @description Generates a heatmap of the distribution of sequences between two populations. In the resulting heatmap, the color of a tile at coordinates (x, y)
#' indicates the number of sequences present in x males, on the horizontal axis, and y females, on the vertical axis.
#'
#' @param data A table of distribution of sequences between sexes obtained with the \code{\link{load_sex_distribution_table}} function.
#'
#' @param title Plot title (default NULL).
#'
#' @param show.significance If TRUE, tiles with significant association with sex are highlighted with the color defined in the significance.color parameter (default TRUE).
#'
#' @param significance.color Color of the border for tiles significantly associated with sex (default "red3").
#'
#' @param significance.threshold P-value threshold to consider a tile significantly associated with sex (default 0.05).
#'
#' @param color.scale.bins A vector of values to use as bins in the color palette (default c(0, 1, 5, 25, 100, 1000)).
#'
#' @param color.scale.colors A vector of two colors used to create the color palette gradient default (c("white", "navyblue")).
#'
#' @return A heatmap stored in a ggplot object.
#'
#' @examples
#' heatmap <- sex_distribution_heatmap(data, title = "Distribution of sequences between sexes",
#'                                     significance.color = "blue", significance.threshold = 0.01,
#'                                     color.scale.bins = c(0, 10, 100, 1000),
#'                                     color.scale.colors = c("white", "red3"))


sex_distribution_heatmap <- function(data, title = NULL,
                                     show.significance = TRUE, significance.color = "red3", significance.threshold = 0.05,
                                     color.scale.bins = c(0, 1, 5, 25, 100, 1000), color.scale.colors = c("white", "navyblue")) {

    # Check that parameters were input correctly
    if (length(color.scale.bins) < 2) {
        stop("Not enough bins in the color scale (minimum 2)")
    }
    if (length(color.scale.colors) != 2) {
        stop("Color scale colors vector should be of length 2")
    }

    # Generate color palette from the color.scale settings
    color_palette <- generate_color_palette(color.scale.bins, color.scale.colors)

    # Associate the corresponding bin to each row of the data frame
    data$data$Bin <- factor(unlist(lapply(data$data$Sequences, function(x) names(color_palette)[tail(which(x >= color.scale.bins), n=1)])), levels = names(color_palette))

    # Check if tile is significant by comparing to significance threshold
    data$data$Signif <- as.factor(data$data$Signif & (data$data$Sequences > 0))

    # Generate the plot
    heatmap <- ggplot2::ggplot(data$data, ggplot2::aes(x = Males, y = Females)) +
        ggplot2::geom_tile(ggplot2::aes(fill = Bin), color = "grey50", size = 0.1) +
        ggplot2::theme_bw() +
        ggplot2::theme(plot.margin = ggplot2::margin(5, 5, 5, 5),
                       panel.border = ggplot2::element_rect(size = 0.5, color = "black"),
                       panel.grid = ggplot2::element_blank(),
                       axis.text = ggplot2::element_text(size = 16, color = "black", face = "bold"),
                       axis.title.x = ggplot2::element_text(size = 18, face = "bold", margin = ggplot2::margin(10, 0, 0, 0)),
                       axis.title.y = ggplot2::element_text(size = 18, face = "bold", margin = ggplot2::margin(0, 10, 0, 0)),
                       legend.margin = ggplot2::margin(0, 0, 10, 0),
                       legend.title = ggplot2::element_text(size = 14, face = "bold"),
                       legend.text = ggplot2::element_text(size = 11),
                       legend.key.height = ggplot2::unit(0.05, "npc"),
                       legend.key.width = ggplot2::unit(0.05, "npc"),
                       legend.key = ggplot2::element_rect(size = 0.5, color = "grey80"),
                       legend.position = "right",
                       legend.text.align = 0) +
        ggplot2::scale_fill_manual(name = "Sequences", breaks = names(color_palette), values = color_palette, labels = names(color_palette), drop = FALSE) +
        ggplot2::scale_x_continuous(name = "Number of males", breaks = seq(0, data$n_males, 5), minor_breaks = seq(0, data$n_males, 1), expand = c(0, 0)) +
        ggplot2::scale_y_continuous(name = "Number of females", breaks = seq(0, data$n_females, 5), minor_breaks = seq(0, data$n_females, 1), expand = c(0, 0))

    # Add highlight to significant tiles if specified
    if (show.significance) {
        heatmap <- heatmap +
            ggplot2::geom_tile(data = data$data, ggplot2::aes(x = Males, y = Females, color = Signif), fill = "NA", size = 0.4) +
            ggplot2::scale_color_manual(name = ggplot2::element_blank(), values = c("TRUE"=significance.color, "FALSE"="NA", color_palette),
                                        breaks = c("TRUE"), labels = c("Signif."))
    }

    # Add title if specified
    if (!is.null(title)) {
        heatmap <- heatmap +
            ggplot2::ggtitle(title) +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 20, face = "bold", margin = margin(0, 0, 10, 0)))
    } else {
        heatmap <- heatmap +
            ggplot2::theme(plot.title = ggplot2::element_blank())
    }
}



#' @title Generate color palette
#'
#' @description Generate a palette vector from the specified bins and colors.
#'
#' @param color.scale.bins A vector of values to use as bins in the color palette.
#'
#' @param color.scale.colors A vector of two colors used to create the color palette gradient.
#'
#' @return A color palette vector.


generate_color_palette <- function(color.scale.bins, color.scale.colors) {

    bin_labels <- c()

    # Create names for all bins except the last. Names are the bin value if the bin has size 1, or an interval if the bin has size > 1
    for (i in 1:(length(color.scale.bins) - 1)) {
        if (color.scale.bins[i] == color.scale.bins[i + 1] - 1) {
            temp <- as.character(color.scale.bins[i])
        } else {
            temp <- paste(as.character(color.scale.bins[i]), as.character(color.scale.bins[i + 1] - 1), sep = "-")
        }
        bin_labels <- c(bin_labels, temp)
    }

    # Last bin is special ("> last bin")
    bin_labels <- c(bin_labels, paste0(">", as.character(color.scale.bins[length(color.scale.bins)])))

    # Generate color palette from the bins and the color values
    get_palette <- colorRampPalette(color.scale.colors)
    palette <- get_palette(length(bin_labels))
    names(palette) <- bin_labels

    return(palette)
}
