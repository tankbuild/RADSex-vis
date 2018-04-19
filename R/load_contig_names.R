#' @title Load a contig names file
#'
#' @description Loads a table of contig names.
#'
#' @param input_file_path Path to a contig names file.
#'
#' @return A vector with contigs as names and contig names as values, or NULL if there is no input file
#'
#' @examples
#' lengths <- load_contig_names("contig_names.tsv")

load_contig_names <- function(input_file_path) {

    if (!is.null(input_file_path)) {

        raw_data <- suppressMessages(readr::read_delim(input_file_path, "\t", escape_double = FALSE, col_names = FALSE,  trim_ws = TRUE))
        data <- raw_data$X2
        names(data) <- raw_data$X1
        return(data)

    } else {

        return(NULL)

    }

}