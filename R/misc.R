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
#' @description Check the integrity and the correctness of a flow cytometry dataset
#' @param fcdataset Flow cytometry dataset to be checked.
#' @return check the integrity and the correctness of a flow cytometry dataset
#'
#' @export
check_IDs <- function(fcdataset) {

  collector <- list()

  for (level_1 in names(fcdataset)) {

    int <- fcdataset[[level_1]]

    for (level_2 in names(int)) {

      collector[[paste(level_1, level_2, sep = "_")]] <- rownames(int[[level_2]])

    }

  }

  result <- all(sapply(collector,
                       FUN = identical, collector[[1]]))

  if (result == TRUE) {
    print("Everything looks fine")
  }else {
    print("Something is not correct with the cell IDs")
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
            print("Something must have went wrong ;)")
          }
        }
        else {
          print(paste0("The parameter '", old , "' is not existing in fcd$expr$", subfolder, "."))
        }

      }
    }
  }
  else {
    print("ERROR: Both vectors have to be the same length.")
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
#' @return df_frequency
#'
#' @export
subset_fcd <- function(fcd, size) {

  random_cells <- sample(rownames(fcd[["expr"]][["orig"]]), size = size)

  condor_filter <- filter_fcd(fcdataset = fcd,
                              cell_ids = random_cells)

  return(condor_filter)
}
