% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_load_and_transform.R
\name{transform_data}
\alias{transform_data}
\title{transform_data}
\usage{
transform_data(keep, transformation, original_data, verbose, cofactor = 5)
}
\arguments{
\item{keep}{Vector of the parameter to keep in the analysis.}

\item{transformation}{transformation to perform.}

\item{original_data}{Original data}

\item{verbose}{Logical, if TRUE the transformation parameters are printed.}

\item{cofactor}{cofactor used for 'arcsinh' transformation, default 5, can be set to 150 for HDFC data.}
}
\value{
transformed flow cytometry dataset
}
\description{
Data transformation, this function runs within the \code{\link{prep_fcd}} wrapper, the logicle tranformation are derived from Cytofkit.
}
\details{
transform_data
}
