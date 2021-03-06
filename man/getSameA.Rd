% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_functions.R
\name{getSameA}
\alias{getSameA}
\title{Generate decision congruence matrix}
\usage{
getSameA(cluster_num = NULL)
}
\arguments{
\item{cluster_num}{number of thought cluster for one option. If 3 then option number is 2 * 3 + 1 = 7.}
}
\value{
A decision congruence matrix.
1 for clusters supporting the same decision; 0 otherwise.
}
\description{
\code{getSameA} generate the decision congruence matrix.
}
