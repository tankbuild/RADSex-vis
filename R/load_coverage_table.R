#' @title Load a coverage table
#'
#' @description Loads a table of coverage obtained with RADSex subset, RADSex signif, or RADSex loci.
#'
#' @param input_file_path Path to a table of coverage.
#'
#' @return A data frame containing the coverage data.
#'
#' @examples
#' data = load_coverage_table(input_file_path)

load_coverage_table <- function(input_file_path) {

    data <- suppressMessages(readr::read_delim(input_file_path, "\t", escape_double = FALSE, col_names = TRUE, trim_ws = TRUE))
    return(data)
}