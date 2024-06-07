#' runPCA
#'
#' @title runPCA
#' @description Run a Principal Component Analysis.
#' @param fcd flow cytometry dataset.
#' @param data_slot name of the data slot to use to calculate the PCA, original data (orig) or harmonized/normalized data (norm).
#' @param seed Seed used for the randomization steps.
#' @param prefix Prefix for the output.
#' @return runPCA
#'
#' @export
runPCA <- function(fcd, data_slot = "orig", seed, prefix = NULL) {

  set.seed(seed)

  if (is.null(prefix)) {

    fcd[["pca"]][[data_slot]] <- prcomp(fcd$expr[[data_slot]])$x

  } else {

    fcd[["pca"]][[paste(prefix, data_slot, sep = "_")]] <- prcomp(fcd$expr[[data_slot]])$x

  }

  return(fcd)

}


#' runUMAP
#'
#' @title runUMAP
#' @description Run a UMAP dimensionality reduction.
#' @param fcd flow cytometry dataset.
#' @param input_type data to use for the calculation of the UMAP, e.g. "expr" or "pca".
#' @param data_slot name of the PCA data slot to use to harmonize. If no prefix was added the, *orig*.
#' @param n_neighbors n_neighbors.
#' @param n_components n_components.
#' @param min_dist min_dist.
#' @param metric metric.
#' @param seed Seed used for the randomization steps.
#' @param prefix Prefix for the output.
#' @param n_threads Number of threads to be used in the UMAP calculation.
#' @param top_PCA Number of PCA used in the UMAP calculation.
#' @param ret_model LOGICAL if the UMAP model should be saved for future projection of the data.
#' @import umap
#' @import Rtsne
#' @return runUMAP
#'
#' @export
runUMAP <- function (fcd,
                     input_type,
                     data_slot,
                     n_neighbors = 15,
                     n_components = 2,
                     min_dist = 0.2,
                     metric = "euclidean",
                     seed,
                     prefix = NULL,
                     n_threads = 32,
                     top_PCA = ncol(fcd[[input_type]][[data_slot]]),
                     ret_model = FALSE) {

  set.seed(seed)

  umap_model <- uwot::umap(X = fcd[[input_type]][[data_slot]][,1:top_PCA],
                           n_neighbors = n_neighbors,
                           n_components = n_components,
                           min_dist = min_dist,
                           metric = metric,
                           n_threads = n_threads,
                           ret_model = ret_model)

  if (ret_model == FALSE) {

    umapMat <- umap_model

  } else {

    umapMat <- umap_model$embedding

  }

  colnames(umapMat) <- c("UMAP1", "UMAP2")

  umap_name <- sub("^_", "" , paste(prefix, input_type, data_slot, sep = "_"))

  if (top_PCA < ncol(fcd[[input_type]][[data_slot]])) {

    suffix <- paste0("top", top_PCA)

    umap_name <- paste(umap_name, suffix, sep = "_")
  }

  fcd[["umap"]][[umap_name]] <- umapMat

  if (ret_model == TRUE) {

    fcd[["extras"]][["umap_model"]] <- umap_model

  }

  return(fcd)

}

#' runDM
#'
#' @title runDM
#' @description Run Diffusion Map dimensionality reduction.
#' @param fcd flow cytometry dataset.
#' @param input_type data to use for the calculation of the UMAP, e.g. "pca" (suggested option).
#' @param data_slot name of the PCA data slot to use to harmonize. If no prefix was added the, *orig*.
#' @param k K used for the analysis.
#' @param seed Seed used for the randomization steps.
#' @param prefix Prefix of the output.
#' @param top_PCA Number of principal components to use for the analysis.
#' @import destiny
#' @import SingleCellExperiment
#' @import slingshot
#' @return runDM
#'
#' @export
runDM <- function(fcd,
                  input_type,
                  data_slot,
                  k = 10,
                  seed,
                  prefix = NULL,
                  top_PCA = ncol(fcd[[input_type]][[data_slot]])) {

  set.seed(91)

  dm <- DiffusionMap(fcd[[input_type]][[data_slot]][,1:top_PCA],
                     vars = NULL,
                     k=k, suppress_dpt = TRUE, verbose=TRUE, n_pcs = NA)

  dm <- cbind(dm$DC1, dm$DC2, dm$DC3)
  colnames(dm) <- c("DC_1", "DC_2", "DC_3")

  dm_name <- sub("^_", "" , paste(prefix, input_type, data_slot, sep = "_"))

  if (top_PCA < ncol(fcd[[input_type]][[data_slot]])) {

    suffix <- paste0("top", top_PCA)

    dm_name <- paste(umap_name, suffix, sep = "_")
  }

  fcd[["diffmap"]][[dm_name]] <- dm

  return(fcd)

}

#' runtSNE
#'
#' @title runtSNE
#' @description Calculate tSNE dimensionality reduction..
#' @param fcd flow cytometry dataset.
#' @param input_type data to use for the calculation, e.g. "pca" (suggested option).
#' @param data_slot name of the PCA data slot to use to harmonize. If no prefix was added the, *orig*.
#' @param perplexity Perplexity used for tSNE calculation (see Rtsne documentation for details).
#' @param seed Seed used for the randomization steps.
#' @param prefix Prefix of the output.
#' @param n_threads Number of threads to be used in the tSNE calculation.
#' @param top_PCA Number of principal components to use for the analysis.
#' @return tSNE cohordinates
#'
#' @export
runtSNE <- function(fcd,
                    input_type,
                    data_slot,
                    perplexity,
                    seed,
                    prefix = NULL,
                    n_threads = 1,
                    top_PCA = ncol(fcd[[input_type]][[data_slot]])) {

  set.seed(seed)

  tSNE_df <- Rtsne(X = fcd[[input_type]][[data_slot]][,1:top_PCA],
                   dims = 2,
                   perplexity = perplexity,
                   check_duplicates = F,
                   verbose = T,
                   num_threads = n_threads,
                   pca = FALSE)

  tSNE_df <- tSNE_df$Y

  colnames(tSNE_df) <- c("tSNE1", "tSNE2")

  rownames(tSNE_df) <- rownames(fcd[[input_type]][[data_slot]])

  tSNE_name <- sub("^_", "" , paste(prefix, input_type, data_slot, sep = "_"))

  if (top_PCA < ncol(fcd[[input_type]][[data_slot]])) {

    suffix <- paste0("top", top_PCA)

    tSNE_name <- paste(tSNE_name, suffix, sep = "_")
  }

  fcd[["tSNE"]][[tSNE_name]] <- tSNE_df

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
