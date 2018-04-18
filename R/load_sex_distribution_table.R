#' @title Load a sex distribution table
#'
#' @description Loads a table of distribution of sequences between sexes obtained with RADSex sex_distribution.
#'
#' @param input_file_path Path to a table of distribution of sequences between sexes.
#'
#' @return A list with the following elements:
#' \item{data}{A data frame of the distribution of sequences between sexes.}
#' \item{n_males}{Number of males in the population}
#' \item{n_females}{Number of females in the population}
#'
#' @examples
#' data <- load_sex_distribution_table("sex_distribution_table.tsv")

load_sex_distribution_table <- function(input_file_path) {

    data <- suppressMessages(readr::read_delim(input_file_path, "\t", escape_double = FALSE, col_names = TRUE, trim_ws = TRUE))
    n_males <- max(data$Males)
    n_females <- max(data$Females)

    return(list(data = data, n_males = n_males, n_females = n_females))
}