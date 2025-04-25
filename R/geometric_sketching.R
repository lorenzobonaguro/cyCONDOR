#' runPCA
#'
#' @title runPCA
#' @description Performs a Principal Component Analysis (PCA) on the expression matrix specified in \code{data_slot}.
#' @param fcd flow cytometry dataset.
#' @param data_slot data slot to use for the calculation, e.g. "orig" or "norm".
#' @param seed A seed is set for reproducibility.
#' @param prefix Optional prefix for the slot name of the output.
#' @param markers Vector of marker names to include or exclude from the calculation according to the discard parameter. See functions \code{\link{used_markers}} and \code{\link{measured_markers}} for the extraction of markers directly from the condor object.
#' @param discard LOGICAL if the markers specified should be included, "F", or excluded, "T", from the calculation. Default = F.
#'
#' @import stats
#' @import dplyr
#'
#' @details
#' The calculation of the PC is based on the function \code{\link[stats]{prcomp}} from the R Stats Package. See the RDocumentation (https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/prcomp) for more details.
#'
#' @return The function returns a fcd with a Principle Components (PC) data frame saved in \code{fcd$pca}. The name of the output consists of the prefix (if given) and the data slot.
#'
#' @export



geometric_sketching <- function(fcd,
                   data_slot = "orig",
                   seed = 91,
                   prefix =  NULL,
                   markers = colnames(fcd$expr[["orig"]]),
                   discard = FALSE, max_cells=10000){
  set.seed(seed)


reticulate::use_condaenv("geosketch",required=T)
reticulate::source_python(here::here("scripts/python/grapheno.py"))


gs<-import("geosketch.gs")
  # see if selected markers are present in condor_object
  for (single in markers){
    if (!single %in% colnames(fcd[["expr"]][["orig"]])){
      stop(paste("ERROR:",single, "not found in present markers."))

    }
  }



  # save rotated values of PCs in "pca"-slot of condor
  PCA_array<-fcd[["pca"]][[paste(sub("^_", "", paste(prefix, data_slot, sep = "_")))]] 



sketch_index = gs(PCA_array, max_cells, replace=False)

  return(fcd)
}
