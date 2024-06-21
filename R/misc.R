#' used_markers
#'
#' @title used_markers
#' @description returns a vector of markers for given input_type and data slot.
#'
#' @param fcd flow cytometry dataset.
#' @param input_type Data for marker extraction, e.g. "pca", "umap", "phenograph", "FlowSOM"
#' @param data_slot Data slot for marker extraction, e.g. "orig" or "norm".
#' @param prefix Prefix of the specific data_slot, if used.
#' @param mute LOGICAL, if output of funtion is wanted (F) or not (T). Default = F.
#'
#' @return marker vector.
#'
#' @export


used_markers <- function(fcd,
                         input_type,
                         prefix = NULL,
                         data_slot,
                         mute = F){
  markers <- fcd$extras$markers[[paste(input_type, sub("^_", "", paste(prefix, data_slot,"markers", sep = "_")), sep = "_")]]
  if (mute == F){

    print(paste("number of measured markers:",length(markers)))
    print(markers)
    return(markers)
  }else{
    return(markers)}
}


#' measured_markers
#'
#' @title measured_markers
#' @description returns a vector of markers for orig condor expr.
#' @param fcd flow cytometry dataset.
#' @return marker vector.
#'
#' @export

measured_markers <- function(fcd){
  markers <- colnames(fcd$expr$orig)
  print(paste("number of measured markers:",length(markers)))
  print(markers)
}


#' filter_fcd
#'
#' @title filter_fcd
#' @description Filter a fcd according to selected cell IDs
#' @param fcdataset Flow cytometry dataset to be filtered.
#' @param cell_ids Row names of the cells to be filtered, should be provided as vector.
#' @return filter the flowframe according to specific cell IDs
#'
#' @export
filter_fcd <- function(fcdataset, cell_ids) {

  new_fcd <- list()

  for (level_1 in names(fcdataset)[names(fcdataset) != "extras"]) {

    int_cont <- fcdataset[[level_1]]

    int_collector <- list()

    for (level_2 in names(int_cont)) {

      int_collector[[level_2]] <- int_cont[[level_2]][cell_ids, ]
    }

    new_fcd[[level_1]] <- int_collector

  }

  if (length(names(fcdataset)[names(fcdataset) == "extras"]) == 1) {

    new_fcd[["extras"]] <- fcdataset[["extras"]]

  }

  return(new_fcd)

}

#' check_IDs
#'
#' @title check_IDs
#' @description Check the integrity and the correctness of a flow cytometry dataset by comparing the cell IDs of all existing slots to the cell IDs of the orig. exprs
#' @param fcd Flow cytometry dataset to be checked.
#' @return check the integrity and the correctness of a flow cytometry dataset
#'
#' @export
#'
check_IDs <- function(fcd){

  collector <- list()
  for (level_1 in names(fcd)) {
    # exclude extra-slot from comparison
    if (!level_1 == "extras"){
      int <- fcd[[level_1]]
      for (level_2 in names(int)) {
        if (!is.null(rownames(int[[level_2]]))){
          # save rownames of levels that are not NULL
          collector[[paste(level_1, level_2, sep = "$")]] <- rownames(int[[level_2]])
        }
        else {

          stop(paste0("The cell IDs of ",level_1, "$",level_2, " are not defined."))

        }
      }
    }}
  result <- all(sapply(collector, FUN = identical, collector[[1]]))
  if (result == TRUE) {
    print("Everything looks fine")
  } else {

    # If something is wrong: print which cell IDs are different
    output <- as.data.frame(sapply(collector, FUN = identical, collector[[1]]))
    colnames(output) <- "comparison_result"
    faulty_IDs_levels <- rownames(output%>%filter(comparison_result==FALSE))
    warning(paste("Something is not correct with the cell IDs in level:",faulty_IDs_levels))

  }
}

#' merge_condor
#'
#' @title merge_condor
#' @description Merges two condor objects.
#' @param data1 Dataset 1 to merge.
#' @param data2 Dataset 2 to merge.
#' @return Condor Object
#'
#' @export
merge_condor <- function(data1, data2) {

  if (identical(colnames(data1$expr$orig), colnames(data2$expr$orig))) {
    paste("Markers are matching")
  } else {
    stop("Markers are not matching")
  }

  if (identical(colnames(data1$anno$cell_anno), colnames(data2$anno$cell_anno))) {
    paste("Metadata are matching")
  } else {
    stop("Metadata are not matching")
  }

  # check if rownames overlap

  if (length(intersect(rownames(data1$expr$orig), rownames(data2$expr$orig))) == 0 ){
    print("Cell IDs are unique")
  } else {
    stop("Cell IDs are not unique")
  }

  new_obj <- list()

  new_obj[["expr"]][["orig"]] <- rbind(data1$expr$orig, data2$expr$orig)
  new_obj[["anno"]][["cell_anno"]] <- rbind(data1$anno$cell_anno, data2$anno$cell_anno)

  class(new_obj) <- "flow_cytometry_dataframe"

  return(new_obj)

}

#' change_param_name
#'
#' @title change_param
#' @description  change parameter names
#' @param fcd flowframe object (condor).
#' @param old_names vector of names that should be changed.
#' @param new_names vector of new names in the same order.
#' @return change_param
#'
#' @export
change_param_name <- function(fcd,
                              old_names,
                              new_names){

  if (length(old_names) == length(new_names)) {

    for (i in 1:length(old_names)) {

      old = old_names[i]
      new = new_names[i]

      for (subfolder in names(fcd$expr)) {

        if (old %in% names(fcd$expr[[subfolder]])) {

          names(fcd$expr[[subfolder]])[which(names(fcd$expr[[subfolder]]) == old)] <- new

          if (new %in% names(fcd$expr[[subfolder]])){

            print(paste0("Changed parameter '", old, "' to '", new, "' in ", subfolder, "."))

          }
          else {
            stop("Something must have went wrong ;)")
          }
        }
        else {
          stop(paste0("The parameter '", old , "' is not existing in fcd$expr$", subfolder, "."))
        }

      }
    }
  }
  else {
    stop("ERROR: Both vectors have to be the same length.")
  }

  return(fcd)
}

#' df_frequency
#'
#' @title df_frequency
#' @description function to get dataframe of frequencies
#' @param classification classification parameters.
#' @param classification_header **optional** chosen header for classification parameters, default = "classification".
#' @param vertical logical, if FALSE frequency on level of classification, default = TRUE.
#' @param groups ** optional ** vector of selected groups to display, default = all.
#' @param condition grouping to be used.
#' @import reshape2
#' @return df_frequency
#'
#' @export
df_frequency <- function(classification,      #condor$clustering$Phenograph_pca_norm_k60$Phenograph
                         classification_header = "classification",
                         condition,            #condor$anno$cell_anno$group
                         vertical = TRUE,
                         groups = NULL) {

  table_occurance <- table(classification, condition)

  i=2

  if (!vertical){

    i=1
  }

  tmp_df_frequencies <- as.data.frame(round(100 * prop.table(table_occurance, i),2))

  result <- dcast(tmp_df_frequencies, classification ~ condition, value = "Freq")

  if (!is.null(groups)){

    keep = append(c("classification"),groups)

    result <- result[,(colnames(result) %in% keep), drop = FALSE]
  }

  colnames(result)[1] <- classification_header

  return(result)
}

#' subset_fcd
#'
#' @title subset_fcd
#' @description Performs a random subset of the fcd
#' @param fcd flowframe object.
#' @param size Numeric: size of the sub-sampling.
#' @param seed A seed is set for reproducibility.
#' @return df_frequency
#'
#' @export
subset_fcd <- function(fcd, size, seed = 91) {

  set.seed(seed)

  random_cells <- sample(rownames(fcd[["expr"]][["orig"]]), size = size)

  condor_filter <- filter_fcd(fcdataset = fcd,
                              cell_ids = random_cells)

  return(condor_filter)
}

#' checkInput
#'
#' @title Helper function within cyCONDOR
#' @description
#' `checkInput()` is a helper function within several cyCONDOR functions for visualization and differential testing.
#' It checks availability if selected 1st level elements (clustering) and 2nd level elements (expr_slot, cluster_slot, cell_anno),
#' and 3rd level variable names in cluster_slot and cell_anno are present.
#' In order to ensure proper functioning, all arguments that are strictly necessary in the parent function and
#' for which no default can be given, should NOT be set to NULL. Please note that this function does not give an error, when group_var, sample_var or pair_var were set to NULL by user in parent function, since these arguments are sometimes optional in the functions.
#' @param fcd flow cytometry dataset
#' @param check_expr_slot logical indicating if expr_slot should be checked
#' @param check_cluster_slot logical indicating if cluster_slot should be checked
#' @param check_cell_anno logical indicating if cell_anno should be checked
#' @param check_reduction logical indicating if reduction_method should be checked
#' @param expr_slot set to NULL (default) if expr_slot is not needed in parent function, set equal to expr_slot if expr_slot string should be checked
#' @param cluster_slot set to NULL (default) if cluster_slot is not needed in parent function, set equal to cluster_slot if cluster_slot string should be checked
#' @param cluster_var set to NULL (default) if cluster_var is not needed in parent function, set equal to cluster_var if cluster_var string should be checked
#' @param reduction_method set to NULL (default) if reduction_method is not needed in parent function, set equal to reduction_method if reduction_method string should be checked
#' @param reduction_slot set to NULL (default) if reduction_slot is not needed in parent function, set equal to reduction_slot if reduction_slot string should be checked
#' @param group_var set to NULL (default) if group_var is not needed in parent function, set equal to group_var if group_var string should be checked
#' @param sample_var set to NULL (default) if sample_var is not needed in parent function, set equal to sample_var if sample_var string should be checked
#' @param pair_var set to NULL (default) if pair_var is not needed in parent function, set equal to pair_var if pair_var string should be checked
#' @returns
#' `checkInput()` returns an stop message if one of the checks fails.
#'
#'@export
checkInput<-function(fcd,
                     check_expr_slot = F,
                     check_cluster_slot = F,
                     check_cell_anno = F,
                     check_reduction = F,
                     expr_slot = NULL,
                     cluster_slot = NULL,
                     cluster_var = NULL,
                     reduction_method = NULL,
                     reduction_slot = NULL,
                     group_var = NULL,
                     sample_var = NULL,
                     pair_var = NULL
){

  ##check if clustering slot is present
  if(check_cluster_slot == T){
    if(is.null(fcd$clustering)){
      stop('no clustering or cell label prediction was performed on this condor object.')
    }
    if(is.null(cluster_slot)){
      stop('cluster_slot needs to be specified to run this function.')
    }
    if(is.logical(cluster_slot)){
      stop('cluster_slot needs to be a string.')
    }
    # if(is.null(fcd$clustering[[cluster_slot]])){
    #   stop('clustering slot "',cluster_slot,'" is not present in clustering')
    # }
    if(!cluster_slot %in% names(fcd$clustering)){
      stop('clustering slot "',cluster_slot,'" is not present in clustering')
    }

    ##check if cluster_var is present
    if(is.null(cluster_var)){
      stop('cluster_var needs to be specified to run this function.')
    }
    if(!is.null(cluster_var)){
      if(!cluster_var %in% colnames(fcd$clustering[[cluster_slot]])){
        stop('cluster_var: column "',cluster_var,'" is not available in specified clustering slot')
      }
    }
  }

  ##check if expr slot is present
  if(check_expr_slot == T){
    if(is.null(expr_slot)){
      stop('expr_slot needs to be specified to run this function.')
    }
    if(is.logical(expr_slot)){
      stop('expr_slot needs to be a string.')
    }
    # if(is.null(fcd$expr[[expr_slot]])){
    #   stop('expr slot "',expr_slot,'" is not present in expr')
    # }
    if(!expr_slot %in% names(fcd$expr)){
      stop('expr slot "',expr_slot,'" is not present in expr')
    }
  }

  ##check if cell_anno slot is present
  if(check_cell_anno == T){
    if(is.null(fcd$anno$cell_anno)){
      stop('"cell_anno" is not present in anno')
    }
  }

  ##check if reduction slot is present
  if(check_reduction == T){
    if(is.null(reduction_method)){
      stop('reduction_method needs to be specified to run this function.')
    }
    if(is.logical(reduction_method)){
      stop('reduction_method needs to be a string.')
    }
    # if(is.null(fcd[[reduction_method]])){
    #   stop('Dimensionality reduction "',reduction_method,'" is not available in fcd object.')
    # }
    if(!reduction_method %in% names(fcd)){
      stop('Dimensionality reduction "',reduction_method,'" is not available in fcd object.')
    }

    ##check if reduction_slot is present
    if(is.null(reduction_slot)){
      stop('reduction_slot needs to be specified to run this function.')
    }
    if(is.logical(reduction_slot)){
      stop('reduction_slot needs to be a string.')
    }
    if(!reduction_slot %in% names(fcd[[reduction_method]])){
      stop('reduction_slot "',reduction_slot,'" is not available for reduction_method ',reduction_method)
    }
  }


  ##check if cell IDs are in right order
  check <- c(check_expr_slot,check_cluster_slot,check_cell_anno,check_reduction)
  if(sum(check) > 1){
    ID_list <- list()
    if(check_expr_slot == T){
      ID_list[["check_expr_slot"]] <- rownames(fcd$expr[[expr_slot]])
    }
    if(check_cluster_slot == T){
      ID_list[["check_cluster_slot"]] <- rownames(fcd$clustering[[cluster_slot]])
    }
    if(check_cell_anno == T){
      ID_list[["check_cell_anno"]] <- rownames(fcd$anno$cell_anno)
    }
    if(check_reduction == T){
      ID_list[["check_reduction"]] <- rownames(fcd[[reduction_method]][[reduction_slot]])
    }

    ##check IDs
    result <- all(sapply(ID_list, FUN = identical, ID_list[[1]]))
    if (result == FALSE) {
      stop("cell IDs of required input data frames do not match.")
    }
  }



  ## check if group_var, sample_var, pair_var are present
  if(!is.null(group_var)){
    if(!group_var %in% colnames(fcd$anno$cell_anno)){
      stop('column "',group_var,'" is not available in cell_anno')
    }}

  if(!is.null(sample_var)){
    if(!sample_var %in% colnames(fcd$anno$cell_anno)){
      stop('sample_var: column "',sample_var,'" is not available in cell_anno')
    }}

  if(!is.null(pair_var)){
    if(!pair_var %in% colnames(fcd$anno$cell_anno)){
      stop('pair_var: column "',pair_var,'" is not available in cell_anno')
    }}

  #print("all good.")
}
