% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/edist.R
\name{edist.matrix}
\alias{edist.matrix}
\title{Lagged Energy Distance for Matrices}
\usage{
\method{edist}{matrix}(
  x,
  sizes,
  group_ids = NULL,
  type = c("1", "2"),
  lag = 0L,
  weights = rep(1, lag + 1),
  a = 1,
  ...
)
}
\arguments{
\item{x}{a \code{matrix}.}

\item{sizes}{numeric vector of sample sizes.}

\item{group_ids}{character vector of names for each grouped time series.}

\item{type}{character determining the type of lagged energy distance.}

\item{lag}{number of lags.}

\item{weights}{optional vector of length \code{lag} + 1 when \code{type} is 2.
Overrides \code{lag} when incompatible.}

\item{a}{index for energy distance. Must be in \eqn{(0, 2]}.}

\item{...}{unused.}
}
\description{
See \link{edist}.
}
