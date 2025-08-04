#' Train Clinical Classifier
#'
#' @title Train a machine learning classifier for the clinical classification of HDFC data.
#' @description This function trains a machine learning classifier for the clinical classification of HDFC data. This function uses the *CytoDx* framework.
#' @param fcd flow cytometry data set.
#' @param input_type data to use for the calculation, e.g. "expr" (suggested option).
#' @param data_slot Name of the data slot to use for the classification, suggested options are "orig" or "norm".
#' @param sample_names Column name of the metadata table containing the samples names.
#' @param classification_variable Vector (same length as number of cells) with the classes to classify (e.g. ctrl/dis).
#' @param family Response type. Must be one of the following: "gaussian","binomial","poisson","multinomial","cox","mgaussian".
#' @param type1 Type of first level prediction. Type of prediction required. Type "link" gives the linear predictors for "binomial", "multinomial", "poisson" or "cox" models; for "gaussian" models it gives the fitted values. Type "response" gives the fitted probabilities for "binomial" or "multinomial", fitted mean for "poisson" and the fitted relative-risk for "cox"; for "gaussian" type "response" is equivalent to type "link".
#' @param type2 Type of second level prediction.
#' @param parallelCore Number of cores to be used.
#' @param reg If elestic net regularization will be used (Default: FALSE).
#' @param seed A seed is set for reproducibility.
#' @details `train_classifier_model()` is a wrapper function around \code{\link[CytoDx]{CytoDx.fit}} implemented in the package *CytoDx*. 
#' The user can specify all the parameters available for the \code{\link[CytoDx]{CytoDx.fit}} functions, arguments description were copied from the documentation of the *CytoDx* package.
#' The function output a *condor* object including the machine learning model in the *extras* slot.
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
                                   seed = 91) {

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
#' @title Predict clinical classification.
#' @description This function uses the model trained with *train_classifier_model* to predict the clinical classification of new samples. This function uses the *CytoDx* framework.
#' @param fcd flow cytometry data set.
#' @param input_type data to use for the calculation, e.g. "expr" (suggested option).
#' @param data_slot Name of the data slot to use for the classification, suggested options are "orig" or "norm".
#' @param sample_names Column name of the metadata table containing the file names.
#' @param model_object flow cytometry data set with the stored classifier model.
#' @param seed A seed is set for reproducibility.
#' @details `predict_classifier()` is a wrapper function around \code{\link[CytoDx]{CytoDx.pred}} implemented in the package *CytoDx*. 
#' The user can specify all the parameters available for the \code{\link[CytoDx]{CytoDx.pred}} functions, arguments description were copied from the documentation of the *CytoDx* package.
#' The function output a *condor* object including the predicted labels in the *extras* slot.
#' @import CytoDx
#' @return predict_classifier
#'
#' @export
predict_classifier <- function(fcd,
                               input_type,
                               data_slot,
                               sample_names = "expfcs_filename",
                               model_object,
                               seed = 91) {

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
