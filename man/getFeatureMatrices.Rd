% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_functions.R
\name{getFeatureMatrices}
\alias{getFeatureMatrices}
\title{Extract features from formatted data}
\usage{
getFeatureMatrices(
  beta,
  format_data,
  deltaM_value,
  option_num,
  print_names = F
)
}
\arguments{
\item{beta}{Variable names to be included in the memory model, e.g., c('sameC','sameA','dist').
Note that you have to include the associated columns in the data.
For example, to fit sameC, you have to include columns with names of 'sameC_1'...'sameC_C'.
To fit sameA, you have to include 'sameA_1',...,'sameA_C'.
In the memory model, we could include sameC (cluster resampling; sameC_value== 1), sameA (decision congruence; sameA_value== 1),
and dist (semantic congruence; dist_value == 1)}

\item{format_data}{This should be a formated data.frame.
sID: participant IDs.
qID: question ID.
cID: option ID being sampled (can take 1 to C).
tNo: thought number in that trial (from 1 to T).
terminate: 0 if continue sampling; 1 if yes; 0 if no (only needed for decision models)}

\item{deltaM_value}{1 minus memory decay; 0 for full decay (marcov property); 1 for no decay;
8 for estimating it as a free parameter}

\item{option_num}{Number of total clusters in each choice (C)}

\item{print_names}{whether to include column and row namaes; default FALSE}
}
\description{
Extract features from formatted data
}
