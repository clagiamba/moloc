% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_moloc.R
\name{process.dataset}
\alias{process.dataset}
\title{process.dataset}
\usage{
process.dataset(d, suffix = NULL, type = NULL, prior_var = 0.15^2)
}
\arguments{
\item{d}{list}

\item{suffix}{"df1" or "df2"}
}
\value{
data.frame with log(abf) or log(bf)
}
\description{
Internal function, process each dataset list for coloc.abf
}
\author{
Chris Wallace
}

