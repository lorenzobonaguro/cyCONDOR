#' learnUMAP
#'
#' @title learnUMAP
#' @description learnUMAP
#' @param fcd XX
#' @param input_type XX
#' @param data_slot XX
#' @param model XX
#' @param n_epochs XX
#' @param prefix XX
#' @param n_threads XX
#' @param seed XX
#' @return learnUMAP
#'
#' @export
learnUMAP <- function(fcd,
                      input_type,
                      data_slot,
                      model,
                      n_epochs = 100,
                      prefix = NULL,
                      n_threads = 32,
                      seed) {

  set.seed(seed)

  umapMat <- uwot::umap_transform(X = fcd[[input_type]][[data_slot]],
                                  model = model,
                                  n_epochs = n_epochs,
                                  n_threads = n_threads)

  colnames(umapMat) <- c("UMAP1", "UMAP2")

  umap_name <- sub("^_", "" , paste(prefix, input_type, data_slot, sep = "_"))

  fcd[["umap"]][[umap_name]] <- umapMat

  return(fcd)

}

#' train_transfer_model
#'
#' @title train_transfer_model
#' @description train_transfer_model
#' @param fcd XX
#' @param input_type XX
#' @param data_slot XX
#' @param label XX
#' @param method XX
#' @param tuneLength XX
#' @param trControl XX
#' @param seed XX
#' @import caret
#' @import randomForest
#' @return train_transfer_model
#'
#' @export
train_transfer_model <- function(fcd,
                                 input_type,
                                 data_slot,
                                 label,
                                 method = "knn",
                                 tuneLength = 5,
                                 trControl = trainControl(method = "cv"),
                                 seed) {

  container <- list()

  set.seed(seed)

  model <- train(x = fcd[[input_type]][[data_slot]],
                 y = factor(label),
                 method = method,
                 tuneLength = tuneLength,
                 trControl = trControl)

  performance <- ggplot(model) +
    geom_errorbar(data = model$results, aes(ymin = Accuracy - AccuracySD, ymax = Accuracy + AccuracySD), width = 0.4) +
    theme_classic(base_size = 15)

  features <- plot(varImp(model))

  container[["lt_model"]] <- model
  container[["performace_plot"]] <- performance
  container[["features_plot"]] <- features

  fcd[["extras"]][["lt_model"]] <- container

  return(fcd)
}

#' predict_labels
#'
#' @title predict_labels
#' @description predict_labels
#' @param fcd XX
#' @param input_type XX
#' @param data_slot XX
#' @param model_object XX
#' @param label XX
#' @param seed XX
#' @return predict_labels
#'
#' @export
predict_labels <- function(fcd,
                           input_type,
                           data_slot,
                           model_object,
                           label = "predicted_labels",
                           seed) {

  set.seed(seed)

  fcd[["clustering"]][[label]] <- data.frame(Description = "predicted",
                                             predicted_label = predict(model_object$extras$lt_model$lt_model,
                                                                       newdata = fcd[[input_type]][[data_slot]]))

  return(fcd)

}
