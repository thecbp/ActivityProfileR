% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bootstrap_U_star_parallel.R
\name{bootstrap_U_star_parallel}
\alias{bootstrap_U_star_parallel}
\title{Parallelized version of the bootstrap_U_star() function}
\usage{
bootstrap_U_star_parallel(data, n_boot = 1000, n_cores = 2)
}
\arguments{
\item{data}{the dataset containing the activity profiles}

\item{n_boot}{number describing how many bootstrap datasets we want to use}

\item{n_cores}{number describing how many cores to use}
}
\value{
a tibble with the results of all of the bootstrapping
}
\description{
Parallelized version of the bootstrap_U_star() function
}
