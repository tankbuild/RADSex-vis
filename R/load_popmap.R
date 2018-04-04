#' @title Load a population map
#'
#' @description Loads a population map from a popmap file.
#'
#' @param input_file_path Path to a population map.
#'
#' @return A vector with individuals as names and sex as values.
#'
#' @examples
#' popmap = load_popmap(input_file_path)

load_popmap <- function(input_file_path) {

    data <- suppressMessages(readr::read_delim(input_file_path, "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE))
    popmap <- data$X2
    names(popmap) <- data$X1
    return(popmap)
}