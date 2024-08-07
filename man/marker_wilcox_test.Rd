% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/differential_testing.R
\name{marker_wilcox_test}
\alias{marker_wilcox_test}
\title{Differential expression testing on cell level}
\usage{
marker_wilcox_test(
  fcd,
  cluster_slot,
  cluster_var,
  group_var,
  group1,
  group2,
  p.adjust.method = "BH",
  marker = NULL,
  min_cells_per_group = 10
)
}
\arguments{
\item{fcd}{flow cytometry data set, that has been subjected to the clustering or cell type label prediction with cyCONDOR}

\item{cluster_slot}{string specifying which clustering slot to use to find variable specified in cluster_var}

\item{cluster_var}{string specifying variable name in cluster_slot that identifies cell population labels to be used (e.g. clusters, metaclusters or predicted labels).}

\item{group_var}{string indicating variable name in cell_anno that should be used as grouping variable. The grouping variable needs to have two or more groupss}

\item{group1}{string indicating group level in group_var that should be used to select cells for group 1}

\item{group2}{string indicating group level in group_var that should be used to select cells for group 2}

\item{p.adjust.method}{p-value adjustment method to use for multiple comparison testing, e.g "BH" (Benjamini-Hochberg, default) or "bonferroni". All available options can be checked in the documentation of the \code{\link[rstatix]{adjust_pvalue}} function from the package \code{rstatix}.}

\item{marker}{(optional) vector of character strings indicating which features in the expression expr_slot should be considered during testing. by default, all features are tested.}

\item{min_cells_per_group}{Minimum number of cells per group required to include a cell population for differential testing.}
}
\value{
A data frame containing test results for each marker and cell population combination - one combination per row.
\itemize{
 \item{cluster} : cell population that was tested
 \item{marker} : feature from expression data that was tested
 \item{group1} : group name (level) of group 1
 \item{group2} : group name (level) of group 2
 \item{n1 / n2} : absolute cell counts in group 1 (n1) and group 2 (n2)
 \item{mean1 / mean2} : mean marker expression in group 1 (mean1) and group 2 (mean2)
 \item{p} : p-value
 \item{p.adj} : adjusted p-value.
 \item{delta_mean} : delta of mean1 and mean2
}
}
\description{
EXPERIEMTAL FEATURE, we advise against using it yet and suggest you have a look at the Vignette "Differential Analysis" for other options.
\code{marker_wilcox_test()} performes a Wilcoxon Rank Sum Test for two groups of cells for each marker and cell population combination.
}
\details{
marker_wilcox_test


The function \code{marker_wilcox_test()} compares two groups of cells for each marker-cell population combination. Expression values will be extracted from expr_slot "orig", containing the transformed data.
In case group_var has more than two levels, the dataset will be subsetted to the two levels specified in group1 and group2.
Wilcoxon Rank Sum Test is performed using the \code{\link[stats]{wilcox.test}} implemented in the \code{stats} package. Afterwards p-value adjustment is performed considering all comparisons that were made.
}
