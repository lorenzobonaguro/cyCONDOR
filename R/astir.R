#' runAstir_celltype
#'
#' @title runAstir_celltype
#' @description Predict cell types using Astir. This package requires the python library `astir` and `reticulate` to work.
#' @param fcd Flow cytometry dataset.
#' @param data_slot Data slot to use for the analysis (e.g. "orig" or "norm").
#' @param analysis_path Full path to the output folder of astir analysis.
#' @param manifest_name Filename of the manifest file, this file must be located in the `analysis_path` folder.
#' @param max_epochs Maximum number of epochs.
#' @param learning_rate Learning Rate.
#' @param initial_epochs Initial Epochs.
#' @import reticulate
#' @return runAstir_celltype
#'
#' @export
runAstir_celltype <- function(fcd,
                              data_slot,
                              analysis_path,
                              manifest_name,
                              max_epochs,
                              learning_rate,
                              initial_epochs) {

  # Save the expression matrix as csv
  write.csv(x = fcd$expr[[data_slot]], file = paste0(analysis_path, "expr.csv"))

  # Import the astir python package
  ast_fun <- reticulate::import("astir")

  # Set the directories for the analysis
  expr_path <- paste0(analysis_path, "expr.csv")

  manifest_path <- paste0(analysis_path, manifest_name)

  ast <- ast_fun$from_csv_yaml(csv_input = expr_path, marker_yaml = manifest_path)

  batch_size = dim(ast$get_type_dataset()$get_exprs_df())[1]/100

  ast$fit_type(max_epochs = as.integer(max_epochs),
               batch_size = as.integer(batch_size),
               learning_rate = learning_rate,
               n_init_epochs = as.integer(initial_epochs))

  print(table(ast$get_celltypes()))

  probabilities <- ast$get_celltype_probabilities()

  diagnostic <- ast$diagnostics_celltype()

  cell_types <- ast$get_celltypes()

  write.csv(x = probabilities, file = paste0(analysis_path, "probabilities.csv"))

  write.csv(x = diagnostic, file = paste0(analysis_path, "diagnostic.csv"))

  write.csv(x = cell_types, file = paste0(analysis_path, "cell_types.csv"))

  # Add the cell type prediction in the condor object

  # Prepare the dataframe
  df <- data.frame(cell_type = cell_types,
                   Description = paste0("Astir_cell_type_", data_slot,
                                        "_Max_Epoc_", max_epochs,
                                        "Learning_Rate_", learning_rate,
                                        "_Initial_Epochs_", initial_epochs))

  fcd[["astir"]][[paste0("Astir_cell_type_", data_slot)]] <- df

  return(fcd)

}

#' runAstir_cellstates
#'
#' @title runAstir_cellstates
#' @description Predict cell states using Astir. This package requires the python library `astir` and `reticulate` to work.
#' @param fcd Flow cytometry dataset.
#' @param data_slot Data slot to use for the analysis (e.g. "orig" or "norm").
#' @param analysis_path Full path to the output folder of astir analysis.
#' @param manifest_name Filename of the manifest file, this file must be located in the `analysis_path` folder.
#' @param max_epochs Maximum number of epochs.
#' @param learning_rate Learning Rate.
#' @param initial_epochs Initial Epochs.
#' @import reticulate
#' @importFrom utils write.csv
#' @return runAstir_cellstates
#'
#' @export
runAstir_cellstates <- function(fcd,
                                data_slot,
                                analysis_path,
                                manifest_name,
                                max_epochs,
                                learning_rate,
                                initial_epochs) {

  # Save the expression matrix as csv
  write.csv(x = fcd$expr[[data_slot]], file = paste0(analysis_path, "expr.csv"))

  # Import the astir python package
  ast_fun <- reticulate::import("astir")

  # Set the directories for the analysis
  expr_path <- paste0(analysis_path, "expr.csv")

  manifest_path <- paste0(analysis_path, manifest_name)

  # Parpare to run astir
  ast <- ast_fun$from_csv_yaml(csv_input = expr_path, marker_yaml = manifest_path)

  batch_size = dim(ast$get_type_dataset()$get_exprs_df())[1]/100

  ast$fit_state(max_epochs = as.integer(max_epochs),
                batch_size = as.integer(batch_size),
                learning_rate = learning_rate,
                n_init_epochs = as.integer(initial_epochs))

  diagnostic <- ast$diagnostics_cellstate()

  cell_states <- ast$get_cellstates()

  write.csv(x = diagnostic, file = paste0(analysis_path, "diagnostic.csv"))

  write.csv(x = cell_states, file = paste0(analysis_path, "cell_states.csv"))

  # Add the cell type prediction in the condor object

  # Prepare the dataframe
  df <- data.frame(cell_states,
                   Description = paste0("Astir_cell_state_", data_slot,
                                        "_Max_Epoc_", max_epochs,
                                        "Learning_Rate_", learning_rate,
                                        "_Initial_Epochs_", initial_epochs))

  fcd[["astir"]][[paste0("Astir_cell_state_", data_slot)]] <- df

  return(fcd)

}
