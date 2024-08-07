% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/batch_correction.R
\name{run_cytonorm}
\alias{run_cytonorm}
\title{run_cytonorm}
\usage{
run_cytonorm(
  fcd,
  batch_var,
  keep_fcs = TRUE,
  output_dir = paste0("./CytoNorm_output_", Sys.Date()),
  files = NULL,
  data_path = NULL,
  anno_table = NULL
)
}
\arguments{
\item{fcd}{flow cytometry dataset}

\item{batch_var}{Column name of batch variable from \code{fcd$anno$cell_anno}.}

\item{keep_fcs}{Boolean whether to keep the normalized FCS files in \code{output_dir}.}

\item{output_dir}{Directory to save normalized FCS files temporary or permanently, if \code{keep_fcs == TRUE}.}

\item{files}{Vector of fcs file names of samples which should be normalized. By default all files contained in the flow cytometry dataset are used.}

\item{data_path}{File path to folder where .fcs files contained in the fcd are stored.}

\item{anno_table}{Path to the annotation table file.}
}
\value{
fcd with a normalized expression data frame saved in \code{fcd$expr$norm}.
}
\description{
Wrapper function around \code{\link[CytoNorm]{CytoNorm.normalize}} from the CytoNorm package.
}
\details{
run_cytonorm


This function assumes that your fcd contains a trained model computed by \code{\link{train_cytonorm}}. The function performs normalization of the samples contained in your fcd. The normalized expression values are added to your fcd and by default FCS files with the normalized values are written to the \code{output_dir}. See [Van Gassen et al., 2019](https://doi.org/10.1002/cyto.a.23904) for more details on CytoNorm.
}
