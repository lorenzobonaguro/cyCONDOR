#' runPseudotime
#'
#' @title Function of calculate Pseudotime from HDFC data analysed with cyCONDOR
#' @description Function of calculate Pseudotime from HDFC data analysed with cyCONDOR. This function uses *slingshot* to calculate trajectories and pesudotime for each cell in the 'condor' object.
#' @param fcd flow cytometry dataset.
#' @param reduction_method Type of dimensionality reduction to calculate the pseudotime.
#' @param reduction_slot Name of the dimensionality reduction slot to be used.
#' @param start.clus (optional) character, indicates the starting cluster(s) from which lineages will be drawn.
#' @param end.clus (optional) character, indicates which cluster(s) will be forced to be leaf nodes in the graph.
#' @param cluster_slot string specifying which clustering slot to use to find variable specified in cluster_var.
#' @param cluster_var string specifying variable name in cluster_slot that identifies cell population labels to be used (e.g. clusters, metaclusters or predicted labels).
#' @param dist.method (optional) character, specifies the method for calculating distances between clusters. Default is "slingshot", see createClusterMST for details.
#' @param use.median logical, whether to use the median (instead of mean) when calculating cluster centroid coordinates.
#' @param omega (optional) numeric or logical, this granularity parameter determines the distance between every real cluster and the artificial cluster, .OMEGA. In practice, this makes omega the maximum allowable distance between two connected clusters. By default, omega = Inf. If omega = TRUE, the maximum edge length will be set to the median edge length of the unsupervised MST times a scaling factor (omega_scale, default = 1.5). This value is provided as a potentially useful rule of thumb for datasets with outlying clusters or multiple, distinct trajectories. See outgroup in createClusterMST.
#' @param omega_scale (optional) numeric, scaling factor to use when omega = TRUE. The maximum edge length will be set to the median edge length of the unsupervised MST times omega_scale (default = 3). See outscale in createClusterMST.
#' @param times numeric, vector of external times associated with either clusters or cells. See defineMSTPaths for details.
#' @param shrink logical or numeric between 0 and 1, determines whether and how much to shrink branching lineages toward their average prior to the split (default = TRUE).
#' @param extend character, how to handle root and leaf clusters of lineages when constructing the initial, piece-wise linear curve. Accepted values are 'y' (default), 'n', and 'pc1'. See 'Details' for more.
#' @param reweight logical, whether to allow cells shared between lineages to be reweighted during curve fitting. If TRUE (default), cells shared between lineages will be iteratively reweighted based on the quantiles of their projection distances to each curve. See 'Details' for more.
#' @param reassign logical, whether to reassign cells to lineages at each iteration. If TRUE (default), cells will be added to a lineage when their projection distance to the curve is less than the median distance for all cells currently assigned to the lineage. Additionally, shared cells will be removed from a lineage if their projection distance to the curve is above the 90th percentile and their weight along the curve is less than 0.1.
#' @param thresh numeric, determines the convergence criterion. Percent change in the total distance from cells to their projections along curves must be less than thresh. Default is 0.001, similar to principal_curve.
#' @param maxit numeric, maximum number of iterations (default = 15), see principal_curve.
#' @param stretch numeric factor by which curves can be extrapolated beyond endpoints. Default is 2, see principal_curve.
#' @param approx_points numeric, whether curves should be approximated by a fixed number of points. If FALSE (or 0), no approximation will be performed and curves will contain as many points as the input data. If numeric, curves will be approximated by this number of points (default = 150 or #cells, whichever is smaller). See 'Details' and principal_curve for more.
#' @param smoother choice of scatter plot smoother. Same as principal_curve, but "lowess" option is replaced with "loess" for additional flexibility.
#' @param shrink.method character denoting how to determine the appropriate amount of shrinkage for a branching lineage. Accepted values are the same as for kernel in density (default is "cosine"), as well as "tricube" and "density". See 'Details' for more.
#' @param allow.breaks logical, determines whether curves that branch very close to the origin should be allowed to have different starting points.
#' @param seed A seed is set for reproducibility.
#' @details `runPseudotime()` is a wrapper function around \code{\link[slingshot]{getLineages}} and \code{\link[slingshot]{getCurves}} implemented in the package *slingshot*.
#' The function first calculated the linages and then the curves of the pseudotime and converts the result to a format compatible with the structure of the 'condor' object. The user can specify all the parameters available for the \code{\link[slingshot]{getLineages}} and  \code{\link[slingshot]{getCurves}} functions, arguments description were copied from the documentation of the *slingshot* package.
#' @return runPseudotime
#' @import slingshot
#' @import DelayedMatrixStats
#'
#' @export
runPseudotime <- function(fcd,
                          reduction_method,
                          reduction_slot,
                          cluster_slot,
                          cluster_var,
                          start.clus = NULL,
                          end.clus = NULL,
                          dist.method = "slingshot",
                          use.median = FALSE,
                          omega = FALSE,
                          omega_scale = 1.5,
                          times = NULL,
                          shrink = TRUE,
                          extend = "y",
                          reweight = TRUE,
                          reassign = TRUE,
                          thresh = 0.001,
                          maxit = 15,
                          stretch = 2,
                          approx_points = NULL,
                          smoother = "smooth.spline",
                          shrink.method = "cosine",
                          allow.breaks = TRUE,
                          seed = 91) {


  ## check input
  checkInput(fcd = fcd,
             check_reduction = T,
             check_cluster_slot = T,
             reduction_method = reduction_method,
             reduction_slot = reduction_slot,
             cluster_slot = cluster_slot,
             cluster_var = cluster_var)

  set.seed(seed)

  # Getting the lineages and convert to data.frame

  print("Slingshot - getLineages")

  lin <- getLineages(data = fcd[[reduction_method]][[reduction_slot]],
                     clusterLabels = fcd$clustering[[cluster_slot]][[cluster_var]],
                     start.clus = start.clus,
                     end.clus = end.clus,
                     dist.method = dist.method,
                     use.median = use.median,
                     omega = omega,
                     omega_scale = omega_scale,
                     times = times)

  lin_df <- slingMST(lin, as.df = TRUE)

  # Getting curves and convert to data.frame

  print("Slingshot - getCurves")

  cur <- getCurves(lin,
                   shrink = shrink,
                   extend = extend,
                   reweight = reweight,
                   reassign = reassign,
                   thresh = thresh,
                   maxit = maxit,
                   stretch = stretch,
                   approx_points = approx_points,
                   smoother = smoother,
                   shrink.method = shrink.method,
                   allow.breaks = allow.breaks)

  cur_df <- slingCurves(cur, as.df = TRUE)

  curve <- as.data.frame(slingPseudotime(cur), row.names = FALSE)

  curve$mean <- rowMeans(curve, na.rm = TRUE)

  rownames(curve) <- rownames(fcd$anno$cell_anno)

  pseudotime_extra <- list(lineages = lin_df, curves = cur_df)

  if (is.null(start.clus)) {

    fcd[["extras"]][[paste("slingshot", reduction_method, reduction_slot, sep = "_")]] <- pseudotime_extra

    fcd[["pseudotime"]][[paste("slingshot", reduction_method, reduction_slot, sep = "_")]] <- curve

  } else {

    fcd[["extras"]][[paste("slingshot", reduction_method, reduction_slot, start.clus, sep = "_")]] <- pseudotime_extra

    fcd[["pseudotime"]][[paste("slingshot", reduction_method, reduction_slot, start.clus, sep = "_")]] <- curve

  }



  return(fcd)

}
