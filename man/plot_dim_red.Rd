% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_visualization.R
\name{plot_dim_red}
\alias{plot_dim_red}
\title{Dimensionality reduction dotplot}
\usage{
plot_dim_red(
  fcd,
  expr_slot = NULL,
  reduction_method,
  reduction_slot,
  cluster_slot = NULL,
  add_pseudotime = FALSE,
  pseudotime_slot,
  param,
  order = FALSE,
  title = "Dimensionality Reduction Plot",
  limX = NULL,
  limY = NULL,
  dot_size = 0.1,
  alpha = 0.2,
  color_discrete = cluster_palette,
  color_gradient = colors,
  remove_guide = FALSE,
  facet_by_variable = FALSE,
  label_clusters = FALSE,
  label_size = 3.5,
  label_color = "black",
  raster = FALSE,
  seed = 91
)
}
\arguments{
\item{fcd}{flow cytometry data set, that has been subjected to dimensionality reduction with cyCONDOR.}

\item{expr_slot}{expr_slot from which to take marker expression values, default is "orig".
Corrected input data should be handled cautiously.}

\item{reduction_method}{string specifying which dimensionality reduction method to use ("umap", "tSNE", "diffmap", "pca").}

\item{reduction_slot}{string specifying reduction name in reduction_method to use for visualization, e.g. "pca_orig".}

\item{cluster_slot}{string specifying which clustering slot to use to find variable specified in cluster_var.}

\item{add_pseudotime}{Logical, if plot should be colored by pseudotime.}

\item{pseudotime_slot}{string specifying pseudotime name to use for visualization.}

\item{param}{parameter to visualize in the plot, this can be either a continuous variable or a categorical one, the function will react differently accordingly.}

\item{order}{logical if you want to order the dots in the plot, by expression for example. This can help to find small populations of positive cells. If set to FALSE, the plotting order of the cells is randomized.}

\item{title}{title of the plot.}

\item{limX}{limits of the x axes (e.g. c(-1, 7)).}

\item{limY}{limits of the y axes (e.g. c(-1, 7)).}

\item{dot_size}{size of the dots.}

\item{alpha}{transparency of the dots.}

\item{color_discrete}{colors for discrete parameters, must be provided as vector of the same length as the number of factors of `param`.}

\item{color_gradient}{colors for continuous parameters.}

\item{remove_guide}{logical, if you want to remove the guide.}

\item{facet_by_variable}{option to facet the plot by a variable, if FALSE the plot is not faceted, if TRUE the plot is faceted by the `param` variable. If any other variable is provided (e.g. "group") the plot will be faceted by this variable.}

\item{label_clusters}{logical: If clusters should be labeled with a text box.}

\item{label_size}{size of the labels.}

\item{label_color}{color of the labels.}

\item{raster}{TRUE or FALSE, if plot should be returned as raster image, this option lowers the quality of the plot but makes it easier to work with images with high number of cells.}

\item{seed}{seed is set for reproducibility.}
}
\value{
plot marker or list of markers
}
\description{
\code{plot_dim_red()} generates a dotplot of the coordinates of any dimensionality reduction performed on a \code{condor} object. The plot can be colored by any variable both numeric (e.g. expression) or categorical (e.g. clustering/metadata).
}
\details{
plot_dim_red
}
