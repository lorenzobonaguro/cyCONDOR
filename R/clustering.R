#' metaclustering
#'
#' @title metaclustering
#' @description Assignment of a metaclusters name.
#' @param fcd flow cytometry dataset.
#' @param clustering Name of the clustering to match for the metaclustering.
#' @param cluster_slot Column containing the original cluster.
#' @param cluster_var Name of the output column.
#' @param metaclusters Vector of the new clusters names, this should be of the same length of the levels of the original clustering.
#' @return metaclustering
#'
#' @export
metaclustering <- function(fcd,
                           clustering,
                           cluster_slot,
                           cluster_var = "metacluster",
                           metaclusters) {

  ### Assign the metacluster
  #### To make it easier when the number of cluster is high

  if (identical(as.character(data.frame(cluster = levels(fcd$clustering[[clustering]][[cluster_slot]]))$cluster), names(metaclusters))) {

    tmp_metaclusters <- data.frame(cluster = levels(fcd$clustering[[clustering]][[cluster_slot]]),
                                   metacluster = c(metaclusters))

    print(tmp_metaclusters)

    fcd$clustering[[clustering]][[cluster_var]] <- factor(fcd$clustering[[clustering]][[cluster_slot]],
                                                       levels = levels(fcd$clustering[[clustering]][[cluster_slot]]),
                                                       labels = tmp_metaclusters$metacluster)

    return(fcd)

  }else {

    print("Names in the metacluster vector do not match the clusters names")

    return(fcd)

  }

}

#' runPhenograph
#'
#' @title runPhenograph
#' @description Run Phenograph based clustering.
#' @param fcd flow cytometry dataset.
#' @param input_type data to use for the calculation, e.g. "expr" or "pca" (suggested: "pca").
#' @param data_slot data slot to use for the calculation, e.g. "orig" or "norm".
#' @param k K value used for clustering.
#' @param seed A seed is set for reproducibility.
#' @param prefix Prefix of the output.
#' @param nPC Number of principal components to use for the analysis.
#' @param markers vector of marker names to include or exclude from the calculation according to the discard parameter. See functions used_markers and measured_markers for the extraction of markers directly from the condor object
#' @param discard LOGICAL if the markers specified should be included, "F", or excluded, "T", from the calculation. Default = F.
#'
#' @import Rphenograph
#' @importFrom igraph membership
#'
#'
#' @return runPhenograph
#'
#' @export


runPhenograph <- function (fcd,
                           input_type,
                           data_slot,
                           k,
                           seed = 91,
                           prefix = NULL,
                           nPC = ncol(fcd$pca[[data_slot]]),
                           markers = colnames(fcd$expr[[data_slot]]),
                           discard = FALSE)
{
  set.seed(seed)

  # see if selected markers are present in condor_object , if input_type = "expr"
  if (input_type == "expr"){
    for (single in markers){
      if (!single %in% colnames(fcd[["expr"]][["orig"]])){
        stop(paste("ERROR:",single, "not found in expr markers."))
      }
    }

    # define markers to use
    if (discard == FALSE){              # (discard == F -> keep specified markers (default = all))

      phe_markers <- markers
    }

    else if (discard == TRUE) {       # (discard == T -> discard specified markers, error if no markers are specified)

      if (length(markers) == length(colnames(fcd$expr[[data_slot]]))){

        stop("ERROR: No markers specified. Specify markers to be removed or set 'discard = F'.")
      }

      else {

        phe_markers <- setdiff(colnames(fcd$expr[[data_slot]]), markers)

      }
    }

    #define fcd subset for calculation
    data1 <- fcd$expr[[data_slot]][, colnames(fcd$expr[[data_slot]]) %in% phe_markers, drop = F]



  }
  if (input_type == "pca"){

    #define fcd subset for calculations and get used markers of PCA analysis

    data1 <- fcd$pca[[data_slot]][,1:nPC]
    phe_markers <- used_markers(fcd,  input_type = "pca", data_slot = data_slot, mute = T)
  }




  Rphenograph_out <- Rphenoannoy::Rphenoannoy(data1, k = k)
  Rphenograph_out <- as.matrix(membership(Rphenograph_out[[2]]))
  Rphenograph_out <- as.data.frame(matrix(ncol = 1,
                                          data = Rphenograph_out,
                                          dimnames = list(rownames(fcd$expr$orig), "Phenograph")))
  Rphenograph_out$Phenograph <- as.factor(Rphenograph_out$Phenograph)
  Rphenograph_out$Description <- paste(input_type, "_", data_slot, "_k", k, sep = "")

  phe_name <- paste("phenograph", sub("^_", "", paste(prefix, input_type, data_slot, "k", k, sep = "_")), sep = "_")

  if (input_type == "pca"){
  if (nPC < ncol(fcd[[input_type]][[data_slot]])) {

    suffix <- paste0("top", nPC)

    phe_name <- paste(phe_name, suffix, sep = "_")
  }
  }

  fcd[["clustering"]][[phe_name]] <- Rphenograph_out

  #save used markers in "extras"-slot
  fcd[["extras"]][["markers"]][[paste(phe_name,"markers", sep = "_")]] <- phe_markers


  return(fcd)
}


#' runFlowSOM
#'
#' @title runFlowSOM
#' @description Run FlowSOM based clustering.
#' @param fcd flow cytometry dataset.
#' @param input_type data to use for the calculation, e.g. "expr" or "pca".
#' @param data_slot data slot to use for the calculation, e.g. "orig" or "norm".
#' @param nClusters Number of final clusters.
#' @param grid_xdim x-axis size of the FlowSOM grid. Default = 10.
#' @param grid_ydim y-axis size of the FlowSOM grid. Default = 10.
#' @param seed A seed is set for reproducibility.
#' @param prefix Prefix of the output.
#' @param ret_model LOGICAL if the model should be saved for future projection of the data. Default = F.
#' @param nPC Number of principal components to use for the analysis. Default = All.
#' @param markers vector of marker names to include or exclude from the calculation according to the discard parameter. See functions used_markers and measured_markers for the extraction of markers directly from the condor object
#' @param discard LOGICAL if the markers specified should be included, "F", or excluded, "T", from the calculation. Default = F.
#'
#' @import FlowSOM
#' @import dplyr
#'
#' @return runFlowSOM
#'
#' @export



runFlowSOM <-  function (fcd,
                         input_type,
                         data_slot,
                         nClusters,
                         grid_xdim= 10,
                         grid_ydim= 10,
                         seed = 91,
                         prefix = NULL,
                         ret_model = FALSE,
                         nPC = ncol(fcd[[input_type]][[data_slot]]),
                         markers = colnames(fcd$expr[[data_slot]]),
                         discard = FALSE)
{
  set.seed(seed)

  # see if selected markers are present in condor_object , if input_type = "expr"
  if (input_type == "expr"){
    for (single in markers){
      if (!single %in% colnames(fcd[["expr"]][["orig"]])){
        stop(paste("ERROR:",single, "not found in expr markers."))
      }
    }

    # define markers to use
    if (discard == FALSE){              # (discard == F -> keep specified markers (default = all))

      som_markers <- markers
    }

    else if (discard == TRUE) {       # (discard == T -> discard specified markers, error if no markers are specified)

      if (length(markers) == length(colnames(fcd$expr[[data_slot]]))){

        stop("ERROR: No markers specified. Specify markers to be removed or set 'discard = F'.")
      }

      else {

        som_markers <- setdiff(colnames(fcd$expr[[data_slot]]), markers)

      }
    }

    #define fcd subset for FlowSOM calculation
    data1 <- fcd$expr[[data_slot]][, colnames(fcd$expr[[data_slot]]) %in% som_markers, drop = F]



  }
  if (input_type == "pca"){

    #define fcd subset for FlowSOM calculations and get used markers of PCA analysis

    data1 <- fcd$pca[[data_slot]][,1:nPC]
    som_markers <- used_markers(fcd,  input_type = "pca", data_slot = data_slot, mute = T)
  }


  som <- FlowSOM::ReadInput(as.matrix(data1),
                            transform = FALSE,
                            scale = FALSE)



  som <- FlowSOM::BuildSOM(som, colsToUse = 1:(ncol(data1)), xdim= grid_xdim, ydim=grid_ydim)  # Question: colsToUse from PCA filtered object?
  som <- FlowSOM::BuildMST(som)
  labels_pre <- som$map$mapping[, 1]


  out <- FlowSOM::metaClustering_consensus(som$map$codes,
                                           k = nClusters,
                                           seed = seed)


  labels <- out[labels_pre]
  labels <- as.factor(labels)

  Description <- paste(input_type, "_", data_slot, "_k", nClusters, sep = "")
  flSOM <- data.frame(FlowSOM = labels, Description = Description)

  rownames(flSOM) <- rownames(fcd$expr$orig)

  SOM_name <- paste("FlowSOM", sub("^_", "", paste(prefix, input_type, data_slot, "k", nClusters, sep = "_")), sep = "_")


  if (input_type == "pca"){
  if (nPC < ncol(fcd[[input_type]][[data_slot]])) {

    suffix <- paste0("top", nPC)
    SOM_name <- paste(SOM_name, suffix, sep = "_")
  }
  }


  if (ret_model == TRUE) {
    fcd[["extras"]][["FlowSOM_model"]] <- som
  }


  fcd[["clustering"]][[SOM_name]] <- flSOM


  #save used markers in "extras"-slot
  fcd[["extras"]][["markers"]][[paste(SOM_name,"markers", sep = "_")]] <- som_markers

  return(fcd)
}
