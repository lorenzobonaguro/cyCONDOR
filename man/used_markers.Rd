% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{used_markers}
\alias{used_markers}
\title{used_markers}
\usage{
used_markers(fcd, input_type, prefix = NULL, data_slot, mute = F)
}
\arguments{
\item{fcd}{flow cytometry dataset.}

\item{input_type}{Data for marker extraction, e.g. "pca", "umap", "phenograph", "FlowSOM"}

\item{prefix}{Optional prefix of the specific data_slot, if used.}

\item{data_slot}{Data slot for marker extraction, e.g. "orig" or "norm".}

\item{mute}{LOGICAL, if output of function is wanted (F) or not (T). Default = F.}
}
\value{
The function returns a vector containing the markers used for the specified method.
}
\description{
returns a vector of markers for given input_type and data slot.
}
\details{
used_markers
}
