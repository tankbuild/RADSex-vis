#' @title Cluster coverage data
#'
#' @description Cluster individuals and sequences in a coverage table.
#'
#' @param data A table of coverage obtained with the \code{\link{load_coverage_table}} function.
#'
#' @param min.coverage Minimum coverage value to consider a sequence present in an individual. Coverage lower than min.coverage will be set to 0.
#'
#' @param max.coverage Maximum sequence coverage allowed in an individual. Coverage higher than max.coverage will be set to max.coverage.
#'
#' @param distance.method Method to use to compute the distance matrix.
#' Possible values: "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski". See \code{\link[stats]{dist}} for details.
#'
#' @param clustering.method Method to use in the clustering.
#' Possible values: "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
#' See \code{\link[stats]{hclust}} for details.
#'
#' @return A list with the following elements:
#' \item{data}{A data frame of coverage with individuals and sequences ordered based on the clustering results}
#' \item{individuals}{Individuals clustering results}
#' \item{sequences}{Sequences clustering results}
#' \item{distribution}{Distribution of coverage values}
#'
#' @examples
#' clustering_data = coverage_clustering(data,
#'                                       min.coverage=0, max.coverage=100,
#'                                       distance.method="binary", clustering.method="complete")


coverage_clustering <- function(data,
                               min.coverage=0, max.coverage=150,
                               distance.method="euclidean", clustering.method="ward.D") {

    # Remove the sequence column which is not useful here
    data <- data[, -2]

    # Count numbers of individuals and sequences
    n_individuals <- dim(data)[2] - 1
    n_sequences <- dim(data)[1]

    # Extract only coverage values for clustering
    coverage <- as.matrix(data[, -1], rownames.force = TRUE)
    rownames(coverage) <- data$id

    # Set lower and upper borders for coverage
    coverage[which(coverage>max.coverage)] <- max.coverage
    coverage[which(coverage<min.coverage)] <- 0
    data[,-1] <- coverage

    # Cluster individuals
    distances <- dist(t(coverage), method = distance.method)
    individual_clusters <- hclust(distances, method = clustering.method)
    data <- data[, c(1, individual_clusters$order + 1)]
    individuals <- individual_clusters$labels[individual_clusters$order]

    # Cluster distances
    distances <- dist(coverage, method = distance.method)
    sequence_clusters <- hclust(distances, method = clustering.method)
    data <- data[sequence_clusters$order,]
    sequences <- sequence_clusters$labels[sequence_clusters$order]

    # Melt data for ggplot
    melted <- suppressMessages(reshape2::melt(data, id.vars = c("id"), variable.name = "individual", value.name = "coverage"))
    melted$id <- factor(melted$id, levels = sequences)
    melted$individual <- factor(as.character(melted$individual), levels = individuals)

    # Compute quantiles for color scale
    distribution <- summary(replace(melted$coverage, which(melted$coverage == 0), NA), na.rm=TRUE)

    output <- list(data = melted, individuals = individual_clusters, sequences = sequence_clusters, distribution = distribution)

    return(output)
}

