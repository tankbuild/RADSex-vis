#' @title Load a mapping results file
#'
#' @description Loads a table of mapping results obtained with the RADSex mapping module.
#'
#' @param input_file_path Path to a mapping results file.
#'
#' @param contig_lengths A vector of contig lengths obtained with the \code{\link{load_contig_lengths} function}.
#'
#' @param contig_names A vector of contig names obtained with the \code{\link{load_contig_names} function}.
#'
#' @param plot.unplaced If TRUE, unplaced contigs will be plotted as a supercontig.
#'
#' @return A list with the following elements:
#' \item{data}{A data frame of the mapping results.}
#' \item{lengths}{A vector with plotted contigs as names and lengths as values}
#'
#' @examples
#' data = load_mapping_results(input_file_path, contig_lengths, contig_names=contig_names, plot.unplaced=FALSE)


load_mapping_results <- function(input_file_path, contig_lengths, contig_names=NULL, plot.unplaced=TRUE) {

    data <- suppressMessages(readr::read_delim(input_file_path, "\t", escape_double = FALSE, trim_ws = TRUE))

    if (!is.null(contig_names)) {  # If contig names are specified, separate the data based on them

        lg_data <- subset(data, data$Contig %in% names(contig_names))
        unplaced_data <- subset(data, !(data$Contig %in% names(contig_names)))

    } else {  # If contig names are not specified, try to separate the data based on contig names starting with LG

        lg_data <- subset(data, substr(data$Contig, 1, 2) == "LG")
        unplaced_data <- subset(data, !(substr(data$Contig, 1, 2) %in% c("LG", "MT")))

    }

    # If lgs were found, sort them and set their color index to 2 (for the plotting later)
    if (dim(lg_data)[1] > 0) {

        lg_data$Contig <- factor(lg_data$Contig, levels = gtools::mixedsort(unique(lg_data$Contig)))
        lg_data$Color <- rep(2, dim(lg_data)[1])
        lengths = c(contig_lengths[which(names(contig_lengths) %in% lg_data$Contig)])

    } else {

        lengths = c()
    }

    if (plot.unplaced & dim(lg_data)[1] != dim(data)[1]) {  # If plot.unplaced is specified and there are unplaced contigs (not just LGs)

        unplaced_data <- unplaced_data[order(unplaced_data$Contig, unplaced_data$Position), ]  # Order by contig first and by position on the contig second

        # Attribute a color index to each unplaced contig, alternating between 0 and 1
        order <- seq(1, length(unique(unplaced_data$Contig)))
        names(order) <- unique(unplaced_data$Contig)
        unplaced_data$Color <- order[unplaced_data$Contig] %% 2

        # Find lengths of unplaced contigs
        unplaced_lengths = contig_lengths[which(names(contig_lengths) %in% names(order))]

        # Transform position on each contig into position on cumulated contig
        temp <- cumsum(unplaced_lengths) - unplaced_lengths[1]
        unplaced_data$Position <- unplaced_data$Position + temp[unplaced_data$Contig]
        unplaced_data$Contig <- "unplaced"

        # Regroup data into one data frame
        data <- rbind(lg_data, unplaced_data)

        lengths = c(lengths, "unplaced"=sum(unplaced_lengths))

    } else {

        # If plot.unplaced was not specified or unplaced contigs were not found, only plot data for LGs
        data <- lg_data

    }

    # Negative log transform p values
    data$P <- -log(data$P, 10)

    # Generate output list
    output <- list(data = data, lengths = lengths)

    return(output)
}