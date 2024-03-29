% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_load_and_transform.R
\name{prepFcsFolderData}
\alias{prepFcsFolderData}
\title{prepFcsFolderData}
\usage{
prepFcsFolderData(LoaderPATH, ceil, useCSV, separator)
}
\arguments{
\item{LoaderPATH}{Path to the .fcs files.}

\item{ceil}{number of cells to subset.}

\item{useCSV}{Logical, if input is .csv and not .fcs.}

\item{separator}{Separator used the flow csv files (if loading from csv).}
}
\value{
load flow cytometry dataset
}
\description{
Load .fcs or .csv files into a dataframe and prepare the condor object.
}
\details{
prepFcsFolderData
}
