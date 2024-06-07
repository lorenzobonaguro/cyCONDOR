#' harmonize_intensities
#'
#' @title harmonize_intensities
#' @description Harmonize the expression values.
#' @param fcd flow cytometry dataset.
#' @param batch vector of column names to use for correcting the data.
#' @param seed Seed used for the randomization steps.
#' @import harmony
#' @return harmonize_intensities
#'
#' @export
harmonize_intensities <- function(fcd, batch, seed) {

  set.seed(seed)

  harmony_param <- HarmonyMatrix(data_mat = as.matrix(fcd$expr$orig),
                                 meta_data = fcd$anno$cell_anno,
                                 vars_use = batch,
                                 do_pca = FALSE)

  fcd[["expr"]][["norm"]] <- as.data.frame(harmony_param)

  return(fcd)

}


#' harmonize_PCA
#'
#' @title harmonize_PCA
#' @description Harmonize the Principal Component Analysis.
#' @param fcd flow cytometry dataset.
#' @param data_slot name of the PCA data slot to use to harmonize. If no prefix was added the, *orig*.
#' @param batch vector of column names to use for correcting the data.
#' @param seed Seed used for the randomization steps.
#' @param prefix Prefix for the output.
#' @return harmonize_PCA
#'
#' @export
harmonize_PCA <- function(fcd, data_slot = "orig", batch, seed, prefix = NULL) {

  set.seed(seed)

  if (is.null(prefix)) {

    fcd[["pca"]][["norm"]] <- HarmonyMatrix(data_mat = as.matrix(fcd$pca[[data_slot]]),
                                            meta_data = fcd$anno$cell_anno,
                                            vars_use = batch,
                                            do_pca = FALSE)

  } else {

    fcd[["pca"]][[paste(prefix, "norm", sep = "_")]] <- HarmonyMatrix(data_mat = as.matrix(fcd$pca[[data_slot]]),
                                                                      meta_data = fcd$anno$cell_anno,
                                                                      vars_use = batch,
                                                                      do_pca = FALSE)

  }

  return(fcd)

}
