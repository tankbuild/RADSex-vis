#' @title Load a contig names file
#'
#' @description Loads a table of contig names.
#'
#' @param input_file_path Path to a contig names file
#'
#' @return A vector with contigs as names and contig names as values.
#'
#' @examples
#' lengths = load_names_file(input_file_path)

load_names_file <- function(input_file_path, group.unplaced=TRUE) {

    raw_data <- suppressMessages(read_delim(input_file_path, "\t", escape_double = FALSE, col_names = FALSE,  trim_ws = TRUE))
    data <- raw_data$X2
    names(data) <- raw_data$X1

    return(data)
}