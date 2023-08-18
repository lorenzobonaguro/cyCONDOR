#' metaclustering
#'
#' @title metaclustering
#' @description Assignment of a metacluster.
#' @param fcd flow cytometry dataset.
#' @param clustering Name of the clustering to match for the metaclustering.
#' @param name_col Column containing the original cluster
#' @param name_out Name of the output column
#' @param metaclusters Vector of the new clusters names, this should be of the same length of the levels of the original clustering.
#' @return metaclustering
#'
#' @export
metaclustering <- function(fcd,
                           clustering,
                           name_col,
                           name_out = "metacluster",
                           metaclusters) {

  ### Assign the metacluster
  #### To make it easier when the number of cluster is high

  if (identical(as.character(data.frame(cluster = levels(fcd$clustering[[clustering]][[name_col]]))$cluster), names(metaclusters))) {

    tmp_metaclusters <- data.frame(cluster = levels(fcd$clustering[[clustering]][[name_col]]),
                                   metacluster = c(metaclusters))

    print(tmp_metaclusters)

    fcd$clustering[[clustering]][[name_out]] <- factor(fcd$clustering[[clustering]][[name_col]],
                                                       levels = levels(fcd$clustering[[clustering]][[name_col]]),
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
#' @description Run Phenograph clustering.
#' @param fcd flow cytometry dataset.
#' @param input_type data to use for the calculation of the UMAP, e.g. "pca" (suggested option).
#' @param data_slot name of the PCA data slot to use to harmonize. If no prefix was added the, *orig*.
#' @param k K value used for clustering.
#' @param seed Seed used for the randomization steps.
#' @param prefix Prefix for the output.
#' @param top_PCA XX
#' @import Rphenograph
#' @importFrom igraph membership
#' @return runPhenograph
#'
#' @export
runPhenograph <- function(fcd,
                          input_type,
                          data_slot,
                          k,
                          seed,
                          prefix = NULL,
                          top_PCA = ncol(fcd[[input_type]][[data_slot]])) {

  set.seed(seed)

  Rphenograph_out <- Rphenoannoy::Rphenoannoy(fcd[[input_type]][[data_slot]][,1:top_PCA], k = k)
  Rphenograph_out <- as.matrix(membership(Rphenograph_out[[2]]))
  Rphenograph_out <- as.data.frame(matrix(ncol=1,data=Rphenograph_out,dimnames=list(rownames(fcd$expr$orig),"Phenograph")))

  Rphenograph_out$Phenograph <- as.factor(Rphenograph_out$Phenograph)
  Rphenograph_out$Description <- paste(input_type, "_",data_slot, "_k", k, sep = "")

  phe_name <- paste("Phenograph", sub("^_", "" , paste(prefix, input_type, data_slot, "k", k, sep = "_")), sep = "_")

  if (top_PCA < ncol(fcd[[input_type]][[data_slot]])) {

    suffix <- paste0("top", top_PCA)

    phe_name <- paste(phe_name, suffix, sep = "_")
  }

  fcd[["clustering"]][[phe_name]] <- Rphenograph_out

  return(fcd)

}

#' runFlowSOM
#'
#' @title runFlowSOM
#' @description Run FlowSOM based clustering.
#' @param fcd flow cytometry dataset.
#' @param input_type data to use for the calculation of the UMAP, e.g. "expr" or "pca".
#' @param data_slot name of the PCA data slot to use to harmonize. If no prefix was added the, *orig*.
#' @param num_clusters number of final clusters.
#' @param seed Seed used for the randomization steps.
#' @param prefix Prefic for the output.
#' @param ret_model XX
#' @param top_PCA XX
#' @return metaclustering
#'
#' @export
runFlowSOM <- function(fcd,
                       input_type,
                       data_slot,
                       num_clusters,
                       seed,
                       prefix = NULL,
                       ret_model = FALSE,
                       top_PCA = ncol(fcd[[input_type]][[data_slot]])) {

  set.seed(seed)

  som <- FlowSOM::ReadInput(as.matrix(fcd[[input_type]][[data_slot]][,1:top_PCA]), transform = FALSE, scale = FALSE)

  som <- FlowSOM::BuildSOM(som, colsToUse = 1:(ncol(fcd[[input_type]][[data_slot]])))

  som <- FlowSOM::BuildMST(som)

  labels_pre <- som$map$mapping[, 1]

  out <- FlowSOM::metaClustering_consensus(som$map$codes, k = num_clusters, seed = seed)

  labels <- out[labels_pre]

  labels <- as.factor(labels)
  Description <- paste(input_type, "_",data_slot, "_k", num_clusters, sep = "")

  flSOM <- data.frame("FlowSOM" = labels, "Description" = Description)

  rownames(flSOM) <- rownames(fcd$expr$orig)

  SOM_name <- paste("FlowSOM", sub("^_", "" , paste(prefix, input_type, data_slot, "k", num_clusters, sep = "_")), sep = "_")

  if (top_PCA < ncol(fcd[[input_type]][[data_slot]])) {

    suffix <- paste0("top", top_PCA)

    SOM_name <- paste(SOM_name, suffix, sep = "_")
  }

  if (ret_model == TRUE) {

    fcd[["extras"]][["FlowSOM_model"]] <- som

  }

  fcd[["clustering"]][[SOM_name]] <- flSOM

  return(fcd)

}


