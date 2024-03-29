% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods-FactorisedExperiment.R,
%   R/methods-ModularExperiment.R, R/methods-ReducedExperiment.R
\name{componentNames<-,FactorisedExperiment-method}
\alias{componentNames<-,FactorisedExperiment-method}
\alias{componentNames<-,ModularExperiment-method}
\alias{moduleNames,ModularExperiment-method}
\alias{moduleNames<-,ModularExperiment-method}
\alias{componentNames,ReducedExperiment-method}
\alias{componentNames<-,ReducedExperiment-method}
\title{Get names of reduced components}
\usage{
\S4method{componentNames}{FactorisedExperiment}(object) <- value

\S4method{componentNames}{ModularExperiment}(object) <- value

\S4method{moduleNames}{ModularExperiment}(object)

\S4method{moduleNames}{ModularExperiment}(object) <- value

\S4method{componentNames}{ReducedExperiment}(object)

\S4method{componentNames}{ReducedExperiment}(object) <- value
}
\arguments{
\item{object}{\link[ReducedExperiment]{ReducedExperiment} object.}
}
\description{
Retrieves the feature names post-dimensionality reduction In the case of
module analysis, these are the names of the gene modules; in the case of
factor analysis, these are the names of the factors.
}
\details{
Component names can be updated with \verb{<-}.
}
