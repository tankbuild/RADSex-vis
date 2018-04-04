#' @title Sex distribution table loader
#'
#' @description Loads a table of distribution of sequences between sexes obtained with RADSex sex_distribution.
#'
#' @param input_file_path Path to a table of distribution of sequences between sexes.
#'
#' @return A list object containing the data, number of males, and number of females.
#'
#' @examples
#' data = load_sex_distribution_table(input_file_path)

load_sex_distribution_table = function(input_file_path) {

    data = suppressMessages(read_delim(input_file_path, "\t", escape_double = FALSE, col_names = TRUE, trim_ws = TRUE))
    n_males = max(data$Males)
    n_females = max(data$Females)

    return(list(data=data, n_males=n_males, n_females=n_females))
}