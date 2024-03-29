% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/samples.R
\name{associate_components}
\alias{associate_components}
\title{Runs linear models for components and sample-level data}
\usage{
associate_components(
  re,
  formula,
  method = "lm",
  scale_reduced = TRUE,
  center_reduced = TRUE,
  type = "II",
  adj_method = "BH",
  ...
)
}
\arguments{
\item{re}{An object inheriting from
\link[ReducedExperiment]{ReducedExperiment}.}

\item{formula}{The model formula to apply. Only the right hand side of the
model need be specified (e.g., "~ x + y"). The left hand side represents the
components themselves. The variables in this formula should be present in the
\code{colData} of \code{re}.}

\item{method}{If "lm", then the \link[stats]{lm} function is used to run
linear models (in tandem with \link[car]{Anova} for running anovas on the
model terms). If "lmer", then linear mixed models are run through
\link[lmerTest]{lmer}.}

\item{scale_reduced}{If TRUE, the reduced data are scaled (to have a standard
deviation of 1) before modelling.}

\item{center_reduced}{If TRUE, the reduced data are centered (to have a mean
of 0) before modelling.}

\item{type}{The type of anova to be applied to the terms of the linear model.}

\item{adj_method}{The method for adjusting for multiple testing. Passed to
the \link[stats]{p.adjust} \code{method} parameter.}

\item{...}{Additional arguments passed to \link[lmerTest]{lmer}, given that
\code{method} is set to "lmer".}
}
\value{
Returns a list with the entry "models" including a list of the
model objects, "anovas" containing the output of anova-based testing,
and "summaries" containing the results of running \code{summary} on the models.
}
\description{
Runs either standard linear or linear mixed models between reduced
components (e.g., factors or modules) and sample-level information.
}
