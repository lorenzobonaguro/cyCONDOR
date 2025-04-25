#' Run Astir cell type prediction
#'
#' @title Subsampling base on geometric sketching
#' @description Subsamples the condor object to a predefined number of cells based on the geosketching algorithm (https://github.com/brianhie/geosketch).
#' Using the PCA representation it subsamples your dataset by preserving the overall structure. Important,
#' if you do not want to loose rare cell populations or do not want skew cell numbers. Requires the python package geosketch and reticulate.
#' @param fcd Flow cytometry dataset.
#' @param pca_slot PCA slot to use for the sketching (e.g. "orig" or "norm").
#' @param n_sub Number of cells to subset (default 50% of all cells).
#' @import reticulate
#' @return a subsampled condor object
#'
#' @export
subsample_geosketch <- function(condor,
                              pca_slot="orig",n_sub=0.5) {


  if( !(pca_slot %in% names(condor$pca)) )
  {
    stop(paste("the slot",pca_slot,"does not exist",sep=" "))
  }

  if (n_sub>nrow(condor$expr$orig)&is.numeric(n_sub)) {
    stop(paste("Subset number",n_sub,"is higher than number of cells",nrow(condor$expr$orig),sep=" "))
  } else if (n_sub<=0) {
    stop(paste("Subset number",n_sub,"needs to be an integer > 0, or a fraction between 0 and 1",sep=" "))
  } else if (!is.integer(n_sub)&!is.numeric(n_sub)&!is.double(n_sub)) {
    stop(paste("Subset number",n_sub,"needs to be an integer > 0, or a fraction between 0 and 1",sep=" "))
  } else if (!(n_sub%%1==0)&!(0<n_sub&n_sub<1)) {
    stop(paste("Subset number",n_sub,"needs to be a fraction between 0 and 1",sep=" "))
  } else if (!(n_sub%%1==0)) {
    n_sub=nrow(condor$expr$orig)*n_sub
  }

  reticulate::use_condaenv("base",required=T)


  gsk<-reticulate::import("geosketch")
  gs<-gsk$gs
  #pd<-reticulate::import("pandas")
  np<-reticulate::import("numpy",convert=F)



  PCA_array=condor$pca[[pca_slot]]
  sketch_index = gs(PCA_array, as.integer(n_sub), replace=F)

  sketch_index_vector<-unlist(sketch_index)

  condor_filter <- filter_fcd(fcd = condor,
                              cell_ids = rownames(condor$expr$orig)[sketch_index_vector])


  return(condor_filter)

}


