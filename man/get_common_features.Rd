% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/features.R
\name{get_common_features}
\alias{get_common_features}
\title{Get common factor features}
\usage{
get_common_features(factor_features)
}
\arguments{
\item{factor_features}{The result of
\link[SummarizedExperiment]{getAlignedFeatures}.}
}
\value{
A data.frame for each factor pair with the numbers and proportions
of the genes in the input that overlap.
}
\description{
Function to count how many genes are aligned with multiple factors.
}
