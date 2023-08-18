#' runPseudotime
#'
#' @title runPseudotime
#' @description Calculate speudotime or flow data
#' @param fcd XX
#' @param dim_red_type XX
#' @param dim_red_name XX
#' @param start.clus XX
#' @param end.clus XX
#' @param clustering XX
#' @param dist.method XX
#' @param use.median XX
#' @param omega XX
#' @param omega_scale XX
#' @param times XX
#' @param shrink XX
#' @param extend XX
#' @param reweight XX
#' @param reassign XX
#' @param thresh XX
#' @param maxit XX
#' @param stretch XX
#' @param approx_points XX
#' @param smoother XX
#' @param shrink.method XX
#' @param allow.breaks XX
#' @param seed XX
#' @import slingshot
#' @import DelayedMatrixStats
#' @return runPseudotime
#'
#' @export
runPseudotime <- function(fcd,
                          dim_red_type,
                          dim_red_name,
                          start.clus = NULL,
                          end.clus = NULL,
                          clustering,
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
                          seed) {

  set.seed(seed)

  # Getting the lineages and convert to data.frame

  print("Slingshot - getLineages")

  lin <- getLineages(data = fcd[[dim_red_type]][[dim_red_name]],
                     clusterLabels = clustering,
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

  fcd[["extras"]][[paste("slingshot", dim_red_type, dim_red_name, sep = "_")]] <- pseudotime_extra

  fcd[["pseudotime"]][[paste("slingshot", dim_red_type, dim_red_name, sep = "_")]] <- curve

  return(fcd)

}
