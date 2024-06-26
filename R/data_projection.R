#' learnUMAP
#'
#' @title learnUMAP
#' @description Projects new samples on a UMAP calculated previously for a reference data set with the same parameters as the new sample. Before executing this function, \code{\link{runUMAP}} needs to be run with \code{ret_model = TRUE} for the reference data set.
#' @param fcd Flow cytometry dataset for which the UMAP coordinates should be predicted.
#' @param input_type Data to use for the calculation of the UMAP, e.g. \code{expr} or \code{pca}. This should be the same which has been used for calculating the UMAP of the reference data set.
#' @param data_slot Name of the \code{input_type} data slot to use e.g. \code{orig}, if no prefix was added. This should be the same which has been used for calculating the UMAP of the reference data set.
#' @param fcd_model Flow cytometry reference data set containing data associated with an existing embedding in \code{fcd_model$extras}.
#' @param nEpochs Number of epochs to use during the optimization of the embedded coordinates. A value between 30 - 100 is a reasonable trade off between speed and thoroughness. By default, this value is set to one third the number of epochs used to build the model.
#' @param prefix Prefix for the name of the dimensionality reduction.
#' @param nThreads Number of threads to use, (except during stochastic gradient descent). By default \code{nThreads = 32}.
#' @param seed A seed is set for reproducibility.
#' @details \code{learnUMAP()} uses \code{\link[uwot]{umap_transform}} to project new samples contained in \code{fcd} on the embedding previously calculated in a reference data set, \code{fcd_model}, using code{\link{runUMAP}}.
#' @return \code{learnUMAP()} returns a \code{fcd} with the predicted UMAP coordinates saved in \code{fcd$umap$expr_orig}, if no \code{prefix} was set.
#'
#' @export
learnUMAP <- function(fcd,
                      input_type,
                      data_slot,
                      fcd_model,
                      nEpochs = 100,
                      prefix = NULL,
                      nThreads = 32,
                      seed = 91) {

  #check if fcd_model contains the model
  if (!"umap_model" %in% names(fcd_model[["extras"]])) {
    stop("Your fcd_model does not contain a 'umap_model' in the 'extras' slot. Rerun 'runUMAP' for your fcd_model with 'ret_model = TRUE'.")
  }

  #check if the column order of fcd and fcd_model
  if(identical(colnames(fcd[[input_type]][[data_slot]]), colnames(fcd_model[[input_type]][[data_slot]]))== FALSE){
    stop("fcd test does not have the same column order as fcd_model")
  }

  set.seed(seed)

  umapMat <- uwot::umap_transform(X = fcd[[input_type]][[data_slot]],
                                  model = fcd_model[["extras"]][["umap_model"]],
                                  n_epochs = nEpochs,
                                  n_threads = nThreads)

  colnames(umapMat) <- c("UMAP1", "UMAP2")

  umap_name <- sub("^_", "" , paste(prefix, input_type, data_slot, sep = "_"))

  fcd[["umap"]][[umap_name]] <- umapMat

  return(fcd)

}

#' train_transfer_model
#'
#' @title train_transfer_model
#' @description Train a machine learning model to transfer cell labels (this function implements the \code{caret} workflow)
#' @param fcd flow cytometry dataset.
#' @param input_type Data to use for the calculation of the UMAP, e.g. \code{expr} or \code{pca}.
#' @param data_slot Name of the \code{input_type} data slot to use e.g. \code{orig}, if no prefix was added.
#' @param cluster_slot string specifying which clustering slot to use to find variable specified in \code{cluster_var}.
#' @param cluster_var string specifying variable name in \code{cluster_slot} that identifies cell population labels to be used (e.g. clusters or metaclusters).
#' @param method A string specifying which classification or regression model to use, by default \code{method = "knn"}. See \code{\link[caret]{train}} for possible values.
#' @param tuneLength An integer denoting the amount of granularity in the tuning parameter grid, default \code{tuneLength = 5}.
#' @param trControl A list of values that define how this function acts, default \code{trControl = caret::trainControl(method = "cv")}. See \code{\link[caret]{trainControl}} and <http://topepo.github.io/caret/using-your-own-model-in-train.html>.
#' @param seed A seed is set for reproducibility.
#' @import caret
#' @import randomForest
#' @details \code{train_transfer_model} uses \code{\link[caret]{train}.
#'
#' @return \code{train_transfer_model} returns a \code{fcd} with the model and associated visualizations saved in \code{fcd$extras$lt_model}.
#'
#' @export
train_transfer_model <- function(fcd,
                                 input_type,
                                 data_slot,
                                 cluster_slot,
                                 cluster_var,
                                 method = "knn",
                                 tuneLength = 5,
                                 trControl = caret::trainControl(method = "cv"),
                                 seed = 91) {

  container <- list()

  set.seed(seed)

  model <- caret::train(x = fcd[[input_type]][[data_slot]],
                        y = factor(fcd[["clustering"]][[cluster_slot]][[cluster_var]]),
                        method = method,
                        tuneLength = tuneLength,
                        trControl = trControl)

  performance <- ggplot(model) +
    geom_errorbar(data = model$results, aes(ymin = Accuracy - AccuracySD, ymax = Accuracy + AccuracySD), width = 0.4) +
    theme_classic(base_size = 15)

  features <- plot(varImp(model))

  container[["lt_model"]] <- model
  container[["performance_plot"]] <- performance
  container[["features_plot"]] <- features

  fcd[["extras"]][["lt_model"]] <- container

  return(fcd)
}

#' predict_labels
#'
#' @title predict_labels
#' @description Uses the model generated with \code{\link{train_transfer_model}} to predict the labels of new samples.
#' @param fcd flow cytometry dataset.
#' @param input_type data to use for the calculation of the UMAP, e.g. "expr" or "pca".
#' @param data_slot name of the PCA data slot to use to harmonize. If no prefix was added the, \code{orig}.
#' @param fcd_model flow cytometry dataset containing Caret model for the label transfer.
#' @param label Label for the output column of the condor object which is saved in the clustering slot of the \code{fcd}.
#' @param seed A seed is set for reproducibility.
#' @return predict_labels
#'
#' @export
predict_labels <- function(fcd,
                           input_type,
                           data_slot,
                           fcd_model,
                           label = "predicted_labels",
                           seed = 91) {

  #check if fcd_model contains the model
  if (!"lt_model" %in% names(fcd_model[["extras"]])) {
    stop("Your fcd_model does not contain 'lt_model' in the 'extras' slot. Run 'train_transfer_model' with your fcd_model first.")
  }

  set.seed(seed)

  fcd[["clustering"]][[label]] <- data.frame(Description = "predicted",
                                             predicted_label = predict(fcd_model$extras$lt_model$lt_model,
                                                                       newdata = fcd[[input_type]][[data_slot]]))

  return(fcd)

}
