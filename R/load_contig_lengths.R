#' @title Load a contig lengths file
#'
#' @description Loads a table of contig lengths.
#'
#' @param input_file_path Path to a contig lengths file.
#'
#' @return A vector with contigs as names and lengths as values.
#'
#' @examples
#' lengths = load_contig_lengths(input_file_path)

load_contig_lengths <- function(input_file_path, group.unplaced=TRUE) {

    raw_data <- suppressMessages(readr::read_delim(input_file_path, "\t", escape_double = FALSE, col_names = FALSE,  trim_ws = TRUE))
    data <- raw_data$X2
    names(data) <- raw_data$X1

    return(data)
}