#' @title Mapping circular plot
#'
#' @description Generates a circular plot of RADSex mapping results in which each sector represents a linkage group and
#' the x-axis represents the position on the linkage group. The y-axis on the first track shows the sex-bias of a sequence, and the second track shows the
#' probability of association with sex of a sequence.
#'
#' @param data A list of mapping results obtained with the \code{\link{load_mapping_results}} function.
#'
#' @param highlights A vector of scaffolds to highlight: c("LG1", "LG7", "scaffold_8") ... (default NULL).
#'
#' @param zoom.highlights If TRUE, highlighted scaffolds will be zoomed on at the top of the plot (default FALSE).
#'
#' @param zoom.ratio Zooming factor for highlighted scaffolds, i.e. a size multiplier (default 2)
#'
#' @param zoom.suffix A suffix to add after the name of a zoomed scaffold (default " (zoom)")
#'
#' @param base.color Background color of a standard sector of the plot (default "white").
#'
#' @param highlight.color Background color of a highlighted sector of the plot (default "grey80").
#'
#' @param point.size Size of a point in the plot (default 0.5)
#'
#' @param color.sex.bias If TRUE, points on the sex-bias track will be colored according to sex.bias.palette (default TRUE).
#'
#' @param sex.bias.palette A vector of three colors defining the sex-bias track palette: female-biased, neutral, male-biased. (default c("firebrick1", "black", "dodgerblue2"))
#'
#' @param color.unmapped If TRUE, unmapped scaffolds will be colored with alternating colors, similar to a manhattan plot (default TRUE).
#'
#' @param unmapped.paltte A named vector of three colors: "0" = alternating color 1, "1" = alternating color2, and "2" = color for mapped scaffolds (default c("0"="dodgerblue3", "1"="goldenrod1", "2"="grey30")).
#'
#' @param signif.threshold Significance threshold for association with sex (default 0.05).
#'
#' @param sector.title.expand A factor that controls the distance between sector titles and sector top axes (default 2).
#'
#' @examples
#'
#' contig_lengths <- load_contig_lengths("contig_lengths.tsv")
#' contig_names <- load_contig_names("contig_names.tsv")
#' data <- load_mapping_results("mapping_results.tsv", contig_lengths, contig_names = contig_names, plot.unplaced = FALSE)
#'
#' mapping_circular_plot(data,
#'                       highlight = NULL, zoom.highlights = FALSE, zoom.ratio = 2, zoom.suffix = " (zoom)",
#'                       base.color = "white", highlight.color = "grey80", point.size = 0.5,
#'                       color.sex.bias = TRUE, sex.bias.palette = c("firebrick1", "black", "dodgerblue2"),
#'                       color.unmapped = TRUE, unmapped.palette = c("0"="dodgerblue3", "1"="goldenrod1", "2"="grey30"),
#'                       signif.threshold = 0.05, sector.titles.expand = 2)


mapping_circular_plot <- function(data,
                                  highlight = NULL, zoom.highlights = FALSE, zoom.ratio = 2, zoom.suffix = " (zoom)",
                                  base.color = "white", highlight.color = "grey80", point.size = 0.5,
                                  color.sex.bias = TRUE, sex.bias.palette = c("firebrick1", "black", "dodgerblue2"),
                                  color.unmapped = TRUE, unmapped.palette = c("0"="dodgerblue3", "1"="goldenrod1", "2"="grey30"),
                                  signif.threshold = 0.05, sector.titles.expand = 1.9) {

    signif.threshold <- -log(signif.threshold / dim(data$data)[1], 10)

    # Get sector information
    n_sectors <- length(data$lengths)
    sectors <- names(data$lengths)
    sector_width <- data$lengths
    if (!is.null(data$names)) {
        sector_names <- data$names
    } else {
        sector_names <- sectors
        names(sector_names) <- sectors
    }

    # Test if lgs specified in "highlight" exist, otherwise remove them from the list and print a warning
    if (!is.null(highlight)) {
        to_remove <- c()
        for (i in 1:length(highlight)) {
            if (!(highlight[i] %in% sectors)) {
                if (!is.null(data$names) & highlight[i] %in% data$names) {
                    highlight[i] <- names(data$names[which(data$names == highlight[i])])
                } else {
                    warning(paste0("Could not find sector \"", highlight[i], "\" given by parameter \"highlight\"."))
                    to_remove <- c(to_remove, i)
                }
            }
        }
        if (!is.null(to_remove)) {highlight <- highlight[-to_remove]}
    }

    # Create the zoomed sector if specified
    if (zoom.highlights & !is.null(highlight)) {

        zoom_data <- data$data[data$data$Contig %in% highlight, , drop = FALSE]  # Extract zoomed sectors data
        zoom_data$Contig <- paste0(zoom_data$Contig, zoom.suffix)  # Edit zoomed sector names with suffix
        data$data <- rbind(data$data, zoom_data)  # Combine base sectors data and zoomed sectors data
        zoom_lengths <- data$lengths[highlight]  # Extract zoomed sector lengths
        names(zoom_lengths) <- paste0(highlight, zoom.suffix)  # Eddit zoomed sector names with suffix
        data$lengths <- c(data$lengths, zoom_lengths)  # Combine base sectors lengths and zoomed sectors lengths
        sectors <- c(sectors, paste0(highlight, zoom.suffix))  # Update sectors list
        n_sectors <- n_sectors + length(highlight)  # Update sectors count
        sector_width <- c(sector_width, zoom_lengths * zoom.ratio)  # Update sector widths
        zoom_names <- paste0(sector_names[highlight], zoom.suffix)
        names(zoom_names) <- paste0(highlight, zoom.suffix)  # Edit zoomed names with suffix
        sector_names <- c(sector_names, zoom_names)
    }

    # Setup sector colors
    bgs <- rep(base.color, n_sectors)
    bgs[which(sectors %in% highlight)] <- highlight.color
    bgs[which(sectors %in% paste0(highlight, zoom.suffix))] <- highlight.color

    # Create sector lengths matrix
    sector_lengths <- matrix(c(rep(0, n_sectors), data$lengths), ncol=2)

    # Setup gaps between sectors
    gaps <- rep(1, n_sectors)
    gaps[length(gaps)] <- 7  # Bigger gap after the last sector to leave some space for track titles
    if (zoom.highlights) { gaps[length(gaps) - length(highlight)] <- 3.5 }  # Add a small gap before the first zoomed sector

    # Assign color to data points (to colorize unmapped scaffold)
    if (color.unmapped) {
        data$data$Color <- unmapped.palette[as.character(data$data$Color)]
    } else {
        data$data$Color <- unmapped.palette["2"]
    }

    data$data$ColorP <- data$data$Color

    if (color.sex.bias) {
        sex_bias_palette <- colorRampPalette(sex.bias.palette)(20)
        data$data$Color[which(data$data$Color == unmapped.palette["2"])] <- sex_bias_palette[(data$data$SexBias[which(data$data$Color == unmapped.palette["2"])] + 1) * 9.5 + 1]
        data$data$ColorP[which(data$data$ColorP == unmapped.palette["2"])] <- unmapped.palette["2"]
    }

    # Calculate angle offset to have zoomed sector on top
    a <- 360 * sum(tail(sector_width, n = length(highlight))) / sum(sector_width) / 2 + 4.75

    circlize::circos.clear()  # Reset circos parameters

    # Setup circos parameters
    circlize::circos.par("track.height" = 0.2,
                         "start.degree" = 90 - a,
                         "gap.degree" = gaps,
                         "cell.padding" = c(0.01, 0.02, 0.01, 0.02),
                         "points.overflow.warning" = FALSE,
                         "track.margin" = c(0, 0.01),
                         "unit.circle.segments" = 100)

    # Initialize circos plot
    circlize::circos.initialize(factors = sectors, xlim = sector_lengths, sector.width = sector_width)

    # Draw the top track of the plot, showing sex-bias
    circlize::circos.track(factors = data$data$Contig,
                           x = data$data$Position,
                           y = data$data$SexBias,
                           ylim = c(-1, 1),
                           bg.col = bgs,
                           panel.fun = function(x, y) {  # panel.fun is the function drawing the track

                               # Get useful sector values
                               sector.index <- circlize::get.cell.meta.data("sector.index")
                               xcenter <- circlize::get.cell.meta.data("xcenter")
                               ymin <- circlize::get.cell.meta.data("ylim")[1]
                               ymax <- circlize::get.cell.meta.data("ylim")[2]
                               xmin <- circlize::get.cell.meta.data("xlim")[1]
                               xmax <- circlize::get.cell.meta.data("xlim")[2]

                               # Add titles to sectors
                               circlize::circos.text(xcenter,
                                                     sector.titles.expand * ymax,
                                                     sector_names[sector.index],
                                                     cex = 1.5,
                                                     facing = "bending.inside",
                                                     niceFacing = TRUE)

                               # Create x axis on top of sectors
                               circlize::circos.axis(h = "top",
                                                     major.at = c(0, xmax/2, xmax),
                                                     labels.cex = 1.2,
                                                     labels.facing = "outside",
                                                     direction="outside",
                                                     labels = convert_to_mb(c(0, xmax/2, xmax)),
                                                     minor.ticks = 4)

                               # Plot the Sex Bias data
                               circlize::circos.points(x, y, cex=point.size, col = data$data$Color[which(data$data$Contig == sector.index)],
                                                       bg = data$data$Color[which(data$data$Contig == sector.index)], pch = 21)

                               # Add a line at 0
                               circlize::circos.segments(xmin, 0, xmax, 0, col = "grey50", lwd = 1)

                               # On the first sector only
                               if (sector.index == sectors[1]) {

                                   # Create y axis
                                   circlize::circos.yaxis(side = "left",
                                                          sector.index = sectors[1],
                                                          labels.cex = 1.2,
                                                          labels.niceFacing = FALSE,
                                                          at = c(1, 0, -1))

                                   # Add y axis labels
                                   # circlize::circos.text(-100,
                                   #                       0.5*(ymax - ymin) + ymin,
                                   #                       "Sex bias",
                                   #                       sector.index = sectors[1],
                                   #                       facing="clockwise",
                                   #                       cex=1.3,
                                   #                       font=2)
                               }
                            }
    )

    # Increase distance between tracks
    circlize::circos.par("track.margin" = c(0, 0.03))

    # Draw the bottom track of the plot, showing -log(p-value association with sex)
    circlize::circos.track(factors = data$data$Contig,
                           x = data$data$Position,
                           y = data$data$P,
                           ylim= round(c(0, max(c(data$data$P, 1.25 * signif.threshold))), 1),
                           bg.col = bgs,
                           panel.fun = function(x, y) { # panel.fun is the function drawing the track

                               # Get useful sector values
                               sector.index <- circlize::get.cell.meta.data("sector.index")
                               xcenter <- circlize::get.cell.meta.data("xcenter")
                               ymin <- circlize::get.cell.meta.data("ylim")[1]
                               ymax <- circlize::get.cell.meta.data("ylim")[2]
                               xmin <- circlize::get.cell.meta.data("xlim")[1]
                               xmax <- circlize::get.cell.meta.data("xlim")[2]

                               # Plot -log(p association with sex)
                               circlize::circos.points(x, y, cex = point.size, col = data$data$ColorP[which(data$data$Contig == sector.index)],
                                                       bg = data$data$ColorP[which(data$data$Contig == sector.index)], pch = 21)

                               # Plot a red line at the significance threshold
                               circlize::circos.segments(xmin, signif.threshold, xmax, signif.threshold, col = "red3", lty = 3, lwd = 1)

                               # Only on first sector
                               if (sector.index == sectors[1]) {

                                   # Create y axis
                                   circlize::circos.yaxis(side = "left",
                                                          sector.index = sectors[1],
                                                          labels.cex = 1.2,
                                                          labels.niceFacing = TRUE,
                                                          at = round(c(0, signif.threshold, max(c(data$data$P, 1.25 * signif.threshold))), 1))

                                   # circlize::circos.text(-0.4*xmax, 0.5*(ymax - ymin) + ymin,
                                   #                       "-log(p)",
                                   #                       sector.index = sectors[1],
                                   #                       facing="clockwise",
                                   #                       cex=1.3,
                                   #                       font=2)
                               }
                           }
    )

    # Decrease distance between tracks
    circlize::circos.par("track.margin" = c(0, 0.01))

    # Add links between base and zoomed sectors
    if (zoom.highlights & !is.null(highlight)) {

        for (i in 1:length(highlight)) {

            circlize::circos.link(highlight[i],
                                  circlize::get.cell.meta.data("cell.xlim", sector.index = highlight[i]),
                                  paste0(highlight[i], zoom.suffix),
                                  circlize::get.cell.meta.data("cell.xlim", sector.index = paste0(highlight[i], zoom.suffix)),
                                  col = rgb(t(col2rgb(highlight.color)), alpha = (255 - 130 * (i-1)/length(highlight)), maxColorValue = 255),
                                  border = rgb(t(col2rgb(highlight.color)), alpha = 255, maxColorValue = 255))
        }
    }
}



# Small function to convert a position value in Mb
convert_to_mb <- function(x) {
    round(x / 10^6, 0)
}

