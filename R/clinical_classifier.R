#' Train Clinical Classifier
#'
#' @title train_classifier_model
#' @description train_classifier_model
#' @param fcd XX
#' @param input_type XX
#' @param data_slot XX
#' @param sample_names XX
#' @param classification_variable XX
#' @param family XX
#' @param type1 XX
#' @param type2 XX
#' @param parallelCore XX
#' @param reg XX
#' @param seed XX
#' @import CytoDx
#' @return train_classifier_model
#'
#' @export
train_classifier_model <- function(fcd,
                                   input_type,
                                   data_slot,
                                   sample_names = "expfcs_filename",
                                   classification_variable,
                                   family = "binomial",
                                   type1 = "response",
                                   type2 = "response",
                                   parallelCore = 1,
                                   reg = FALSE,
                                   seed) {

  set.seed(seed)

  # Perfroms rank transformation
  x_train <- suppressMessages(CytoDx::pRank(x=fcd[[input_type]][[data_slot]],xSample=fcd[["anno"]][["cell_anno"]][[sample_names]]))

  # Convert data frame into matrix. Here we included the 2-way interactions.
  x_train <- model.matrix(~.*.,x_train)

  # Build predictive model using the CytoDx.fit function
  fit <- CytoDx::CytoDx.fit(x=x_train,
                            y=classification_variable,
                            xSample=fcd$anno$cell_anno$expfcs_filename,
                            family = family,
                            type1 = type1,
                            type2 = type2,
                            parallelCore = parallelCore,
                            reg = FALSE)

  fcd[["extras"]][["classifier_model"]] <- fit

  return(fcd)

}


#' Predict Clinical Classifier
#'
#' @title predict_classifier
#' @description predict_classifier
#' @param fcd XX
#' @param input_type XX
#' @param data_slot XX
#' @param sample_names XX
#' @param model_object XX
#' @param seed XX
#' @import CytoDx
#' @return predict_classifier
#'
#' @export
predict_classifier <- function(fcd,
                               input_type,
                               data_slot,
                               sample_names = "expfcs_filename",
                               model_object,
                               seed) {

  set.seed(seed)

  # Perfroms rank transformation
  x_test <- suppressMessages(CytoDx::pRank(x=fcd[[input_type]][[data_slot]], xSample=fcd[["anno"]][["cell_anno"]][[sample_names]]))

  # Convert data frame into matrix. Here we included the 2-way interactions.
  x_test <- model.matrix(~.*.,x_test)

  # Predict AML using CytoDx.ped function
  pred <- CytoDx::CytoDx.pred(model_object,
                              xNew=x_test,
                              xSampleNew=fcd[["anno"]][["cell_anno"]][[sample_names]])

  fcd[["extras"]][["classifier_prediction"]] <- pred

  return(fcd)

}
