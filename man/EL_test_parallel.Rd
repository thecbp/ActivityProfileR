% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EL_test_parallel.R
\name{EL_test_parallel}
\alias{EL_test_parallel}
\title{Paralllelized version of the EL_test function}
\usage{
EL_test_parallel(
  data,
  activity_col,
  group_col,
  grid_length = 100,
  n_boot = 1000,
  quantiles = c(0.05, 0.95),
  alpha = 0.95,
  n_cores = 2
)
}
\arguments{
\item{data}{dataset containing the activity profiles}

\item{grid_length}{number describing how long the activity grid is}

\item{n_boot}{number describing how many bootstrap datasets to make}

\item{quantiles}{two percentiles describing what percentiles of activity points to keep in grid}

\item{alpha}{number describing what quantile to check for the boostrapped test statistics}

\item{n_cores}{number describing how many cores to use for the parallelization}
}
\value{
a list contaning the results of the functional ANOVA test
}
\description{
Paralllelized version of the EL_test function
}
