#' runPCA
#'
#' @title runPCA
#' @description Run a Principal Component Analysis.
#' @param fcd flow cytometry dataset.
#' @param data_slot data slot to use for the calculation, e.g. "orig" or "norm".
#' @param seed A seed is set for reproducibility.
#' @param prefix Prefix of the output.
#' @param markers vector of marker names to include or exclude from the calculation according to the discard parameter. See functions used_markers and measured_markers for the extraction of markers directly from the condor object
#' @param discard LOGICAL if the markers specified should be included, "F", or excluded, "T", from the calculation. Default = F.
#'
#' @import stats
#' @import dplyr
#'
#' @return runPCA
#'
#' @export




runPCA <- function(fcd,
                   data_slot = "orig",
                   seed = 91,
                   prefix =  NULL,
                   markers = colnames(fcd$expr[["orig"]]),
                   discard = FALSE){
  set.seed(seed)

  # see if selected markers are present in condor_object
  for (single in markers){
    if (!single %in% colnames(fcd[["expr"]][["orig"]])){
      stop(paste("ERROR:",single, "not found in present markers."))

    }
  }

  # calculate PC
  if (discard == FALSE){              # (discard == F -> keep markers (default = all))

    tmp <- prcomp(fcd$expr[[data_slot]][, colnames(fcd$expr[["orig"]]) %in% markers, drop = F])

  } else if (discard == TRUE) {       # (discard == T -> discard markers, error if markers are not specified)

    if (length(markers) == length(colnames(fcd$expr[["orig"]]))){     # error code if no markers for removal are specified

      stop("No markers specified. Specify markers to be removed or set 'discard = F'.")


    } else {

      tmp <- prcomp(fcd$expr[[data_slot]][, !colnames(fcd$expr[["orig"]]) %in% markers, drop = F])

    }
  }

  # save rotated values of PCs in "pca"-slot of condor
  fcd[["pca"]][[paste(sub("^_", "", paste(prefix, data_slot, sep = "_")))]] <- tmp$x

  #save used markeres in "extras$markers"-slot
  fcd[["extras"]][["markers"]][[paste("pca", sub("^_", "", paste(prefix, data_slot,"markers", sep = "_")), sep = "_")]] <- dimnames(tmp$rotation)[[1]]



  return(fcd)
}



#' runUMAP
#'
#' @title runUMAP
#' @description Run a UMAP dimensionality reduction.
#'
#' @param fcd flow cytometry dataset.
#' @param input_type data to use for the calculation of the UMAP, e.g. "expr" or "pca".
#' @param data_slot data slot to use for the calculation of the UMAP, e.g. "orig" or "norm".
#' @param nNeighbors Number of neighbors for UMAP calculation. Default = 15.
#' @param nComponents Number of components for UMAP calculation. Default = 2.
#' @param min_dist Min_dist for UMAP calculation. Default = 0.2.
#' @param metric Metric for UMAP calculation. Default = "euclidian". **(other options?)**
#' @param seed A seed is set for reproducibility.
#' @param prefix Prefix of the output.
#' @param nThreads Number of threads to be used in the UMAP calculation. Default = 32.
#' @param nPC Number of PCs used in the UMAP calculation. Default = All.
#' @param ret_model LOGICAL if the UMAP model should be saved for future projection of the data.
#' @param markers vector of marker names to include or exclude from UMAP calculation according to the discard parameter. See functions used_markers and measured_markers for the extraction of markers directly from the condor object
#' @param discard LOGICAL to decide if the markers specified should be included, "F", or excluded, "T", from the UMAP calculation. Default = F.
#'
#' @import umap
#' @import Rtsne **(needed?)**
#'
#' @return runUMAP
#'
#' @export

runUMAP <- function(fcd,
                         input_type,
                         data_slot,
                         nNeighbors = 15,
                         nComponents = 2,
                         min_dist = 0.2,
                         metric = "euclidean",
                         seed = 91,
                         prefix = NULL,
                         nThreads = 32,
                         nPC = ncol(fcd[[input_type]][[data_slot]]),
                         ret_model = FALSE,
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

      UMAP_markers <- markers
    }

    else if (discard == TRUE) {       # (discard == T -> discard specified markers, error if no markers are specified)

      if (length(markers) == length(colnames(fcd$expr[[data_slot]]))){

        stop("ERROR: No markers specified. Specify markers to be removed or set 'discard = F'.")

      }

      else {

        UMAP_markers <- setdiff(colnames(fcd$expr[[data_slot]]), markers)

      }
    }

    #define fcd subset for UMAP calculation
    data1 <- fcd$expr[[data_slot]][, colnames(fcd$expr[[data_slot]]) %in% UMAP_markers, drop = F]



  }

  if (input_type == "pca") {

    #define fcd subset for UMAP calculations and get used markers of PCA analysis
    data1 <- fcd$pca[[data_slot]][,1:nPC]
    UMAP_markers <- used_markers(fcd,  input_type = "pca", data_slot = data_slot, mute = T)
  }
  # calculate UMAP (changed reference to x)
  umap_model <- uwot::umap(X = data1,
                           n_neighbors = nNeighbors,
                           n_components = nComponents,
                           min_dist = min_dist,
                           metric = metric,
                           n_threads = nThreads,
                           ret_model = ret_model)



  # original code no changes
  if (ret_model == FALSE) {
    umapMat <- umap_model
  }
  else {
    umapMat <- umap_model$embedding
  }

  colnames(umapMat) <- c("UMAP1", "UMAP2")


  #name the object
  umap_name <- sub("^_", "", paste(prefix, input_type, data_slot, sep = "_"))

  if (nPC < ncol(fcd[[input_type]][[data_slot]])) {
    suffix <- paste0("top", nPC)
    umap_name <- paste(umap_name, suffix, sep = "_")
  }

  #save the UMAP and model
  fcd[["umap"]][[umap_name]] <- umapMat

  if (ret_model == TRUE) {
    fcd[["extras"]][["umap_model"]] <- umap_model
  }

  #save used markers in "extras"-slot
  fcd[["extras"]][["markers"]][[paste("umap", paste(umap_name,"markers", sep = "_"), sep = "_")]] <- UMAP_markers


  return(fcd)
}

#' runDM
#'
#' @title runDM
#' @description Run Diffusion Map dimensionality reduction.
#' @param fcd flow cytometry dataset.
#' @param input_type data to use for the calculation, e.g. "expr" or "pca".
#' @param data_slot data slot to use for the calculation, e.g. "orig" or "norm".
#' @param k K used for the analysis. Default = 10.
#' @param seed A seed is set for reproducibility.
#' @param prefix Prefix of the output.
#' @param nPC Number of principal components to use for the analysis. Default = All.
#' @param markers vector of marker names to include or exclude from DM calculation according to the discard parameter. See functions used_markers and measured_markers for the extraction of markers directly from the condor object
#' @param discard LOGICAL if the markers specified should be included, "False", or excluded, "True", from the DM calculation. Default = F.
#'
#' @import destiny
#' @import SingleCellExperiment
#' @import slingshot
#' @import dplyr
#'
#' @return runDM
#'
#' @export

runDM <- function (fcd,
                   input_type,
                   data_slot,
                   k = 10,
                   seed = 91,
                   prefix = NULL,
                   nPC = ncol(fcd[[input_type]][[data_slot]]),
                   markers = colnames(fcd$expr[[data_slot]]),
                   discard = FALSE)
{
  set.seed(91)

  # see if selected markers are present in condor_object , if input_type = "expr"
  if (input_type == "expr"){
    for (single in markers){
      if (!single %in% colnames(fcd$expr[[data_slot]])){
        stop(paste("ERROR:",single, "not found in expr markers."))

      }
    }

    # define markers to use
    if (discard == FALSE){              # (discard == F -> keep specified markers (default = all))

      dm_markers <- markers
    }

    else if (discard == TRUE) {       # (discard == T -> discard specified markers, error if no markers are specified)

      if (length(markers) == length(colnames(fcd$expr[[data_slot]]))){

        stop("ERROR: No markers specified. Specify markers to be removed or set 'discard = F'.")

      }

      else {

        dm_markers <- setdiff(colnames(fcd$expr[[data_slot]]), markers)

      }
    }

    #define fcd subset for DM calculation
    data1 <- fcd$expr[[data_slot]][, colnames(fcd$expr[[data_slot]]) %in% dm_markers, drop = F]



  }

  if (input_type == "pca") {

    #define fcd subset for DM calculations and get used markers of PCA analysis
    data1 <- fcd$pca[[data_slot]][,1:nPC]
    DM_markers <- used_markers(fcd,  input_type = "pca", data_slot = data_slot, mute = T)
  }

  # calculate DM (changes to fcd (now data1))

  dm <- DiffusionMap(data1,
                     vars = NULL,
                     k = k,
                     suppress_dpt = TRUE,
                     verbose = TRUE,
                     n_pcs = NA)

  dm <- cbind(dm$DC1, dm$DC2, dm$DC3)

  colnames(dm) <- c("DC_1", "DC_2", "DC_3")

  #name the object
  dm_name <- sub("^_", "", paste(prefix, input_type, data_slot,sep = "_"))

  if (nPC < ncol(fcd[[input_type]][[data_slot]])) {

    suffix <- paste0("top", nPC)

    dm_name <- paste(dm_name, suffix, sep = "_")

  }

  #save the DM
  fcd[["diffmap"]][[dm_name]] <- dm

  #save used markers in "extras"-slot
  fcd[["extras"]][["markers"]][[paste("diffmap", paste(dm_name,"markers", sep = "_"), sep = "_")]] <- dm_markers

  return(fcd)
}

#' runtSNE
#'
#' @title runtSNE
#' @description Calculate tSNE dimensionality reduction.
#' @param fcd flow cytometry dataset.
#' @param input_type data to use for the calculation, e.g. "expr" or "pca" (suggested: "pca").
#' @param data_slot data slot to use for the calculation, e.g. "orig" or "norm".
#' @param perplexity Perplexity used for tSNE calculation (see Rtsne documentation for details).
#' @param seed A seed is set for reproducibility.
#' @param prefix Prefix of the output.
#' @param nThreads Number of threads to be used in the tSNE calculation.
#' @param nPC Number of principal components to use for the analysis.
#' @param markers vector of marker names to include or exclude from UMAP calculation according to the discard parameter. See functions used_markers and measured_markers for the extraction of markers directly from the condor object
#' @param discard Boolean to decide if the markers specified should be included, "F", or excluded, "T", from the UMAP calculation. Default = F.
#'
#' @return tSNE cohordinates
#'
#' @export

runtSNE <- function (fcd,
                         input_type, # expr o. pca
                         data_slot,  # orig, norm or "prefix"-orig/nrom
                         perplexity,
                         seed = 91,
                         prefix = NULL, # new prefix for tSNE DimRed
                         nThreads = 1,
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

      tSNE_markers <- markers
    }

    else if (discard == TRUE) {       # (discard == T -> discard specified markers, error if no markers are specified)

      if (length(markers) == length(colnames(fcd$expr[[data_slot]]))){

        stop("ERROR: No markers specified. Specify markers to be removed or set 'discard = F'.")

      }

      else {

        tSNE_markers <- setdiff(colnames(fcd$expr[[data_slot]]), markers)

      }
    }

    #define fcd subset for tSNE calculation
    data1 <- fcd$expr[[data_slot]][, colnames(fcd$expr[[data_slot]]) %in% tSNE_markers, drop = F]



  }
  if (input_type == "pca"){

    #define fcd subset for tSNE calculations and get used markers of PCA analysis

    data1 <- fcd$pca[[data_slot]][,1:nPC]
    tSNE_markers <- used_markers(fcd,  input_type = "pca", data_slot = data_slot, mute = T)
  }

  tSNE_df <- Rtsne(X = data1,
                   dims = 2,
                   perplexity = perplexity,
                   check_duplicates = F,
                   verbose = T,
                   num_threads = nThreads,
                   pca = FALSE)

  tSNE_df <- tSNE_df$Y
  colnames(tSNE_df) <- c("tSNE1", "tSNE2")
  rownames(tSNE_df) <- rownames(fcd[[input_type]][[data_slot]])

  tSNE_name <- sub("^_", "", paste(prefix, input_type, data_slot, sep = "_"))

  if (nPC < ncol(fcd[[input_type]][[data_slot]])) {

    suffix <- paste0("top", nPC)

    tSNE_name <- paste(tSNE_name, suffix, sep = "_")

  }

  fcd[["tSNE"]][[tSNE_name]] <- tSNE_df

  #save used markers in "extras"-slot
  fcd[["extras"]][["markers"]][[paste("tSNE", paste(tSNE_name,"markers", sep = "_"), sep = "_")]] <- tSNE_markers


  return(fcd)
}

#' runPCA_pseudobulk
#'
#' @title runPCA_pseudobulk
#' @description run a Pricipal component analysis on pseudobulk samples
#' @param fcd flow cytometry dataset.
#' @return runPCA_pseudobulk
#'
#' @export
runPCA_pseudobulk <- function(fcd) {

  data <- fcd[["expr"]][["orig"]]

  index <- fcd[["anno"]][["cell_anno"]]

  container <- list()

  for (sample in unique(index$expfcs_filename)) {

    data_f <- data[rownames(index[index$expfcs_filename == sample,]), ]

    container[[sample]] <- colMeans(data_f)
  }

  data_pca <- do.call(rbind, container)

  pca <- as.data.frame(prcomp(data_pca)$x)

  pca$expfcs_filename <- rownames(pca)

  index <- index[!duplicated(index$expfcs_filename),]

  pca <- merge(pca, index, by = "expfcs_filename")

  output <- list(pca = pca, data = data_pca)

  return(output)

}
