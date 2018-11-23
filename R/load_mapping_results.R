#' @title Load a mapping results file
#'
#' @description Loads a table of mapping results obtained with the RADSex mapping module.
#'
#' @param input_file_path Path to a mapping results file from RADSex map.
#'
#' @param contig_lengths_file_path Path to a contig lengths file from RADSex map.
#'
#' @param contig_names_file_path Path to a file containing the chromosomes names.
#'
#' @param plot.unplaced If TRUE, unplaced scaffolds will be grouped together and plotted as "Unplaced" (default: TRUE).
#'
#' @return A list with the following elements:
#' \item{data}{A data frame of the mapping results.}
#' \item{lengths}{The contig lengths.}
#' \item{names}{A vector with reference contigs names as names and corresponding contig names as values}
#'
#' @examples
#' data <- load_mapping_results("mapping_results.tsv", contig_lengths_file_path = 'contig_lengths.tsv',
#'                              contig_names_file_path = 'contig_names.tsv')


load_mapping_results <- function(input_file_path, contig_lengths_file_path, contig_names_file_path = NULL, plot.unplaced = TRUE) {

    output <- list()

    print(" - Loading contig names file")
    output$names <- load_contig_names(contig_names_file_path)

    print(" - Loading contig lengths file")
    output$lengths <- load_contig_lengths(contig_lengths_file_path, contig_names = output$names, plot.unplaced = plot.unplaced)

    data <- suppressMessages(readr::read_delim(input_file_path, "\t", escape_double = FALSE, trim_ws = TRUE))
    data_lg <- subset(data, data$Contig %in% names(output$lengths$lg))
    data_unplaced <- subset(data, data$Contig %in% names(output$lengths$unplaced))

    # If lgs were found, sort them and set their color index to 2 (for the plotting later)
    if (dim(data_lg)[1] > 0) {

        data_lg$Color <- rep(2, dim(data_lg)[1])

    }

    if (plot.unplaced & dim(data_unplaced)[1] > 0) {  # If plot.unplaced is specified and there are unplaced contigs (not just LGs)

        # Order unplaced contigs data by contig length and then by position on the contig
        data_unplaced <- data_unplaced[order(match(data_unplaced$Contig, names(output$lengths$unplaced)), data_unplaced$Position), ]

        # Attribute a color index to each unplaced contig, alternating between 0 and 1
        order <- seq(1, length(unique(data_unplaced$Contig)))
        names(order) <- unique(data_unplaced$Contig)
        data_unplaced$Color <- order[data_unplaced$Contig] %% 2

        # Transform position on each contig into position on cumulated contig
        temp <- cumsum(output$lengths$unplaced) - output$lengths$unplaced[1]
        data_unplaced$Original_position <- data_unplaced$Position
        data_unplaced$Position <- data_unplaced$Position + temp[data_unplaced$Contig]
        data_unplaced$Contig_id <- data_unplaced$Contig
        data_unplaced$Contig <- "Unplaced"

        # Regroup data into one data frame
        data_lg$Contig_id <- data_lg$Contig
        data_lg$Original_position <- data_lg$Position
        data <- rbind(data_lg, data_unplaced)
        data$Contig <- factor(data$Contig, levels = c(names(output$lengths$lg), "Unplaced"))

    } else {

        # If plot.unplaced was not specified or unplaced contigs were not found, only plot data for LGs
        data <- data_lg
        data$Contig <- factor(data$Contig, levels = names(output$lengths$lg))
        data$Contig_id <- data_lg$Contig
        data$Original_position <- data_lg$Position

    }

    data$P <- -log(data$P, 10)
    output$data <- data

    return(output)
}
