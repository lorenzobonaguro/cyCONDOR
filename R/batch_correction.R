#' harmonize_intensities
#'
#' @title harmonize_intensities
#' @description Harmonize the expression values.
#' @param fcd flow cytometry dataset.
#' @param batch_var vector of column names from \code{fcd$anno$cell_anno} to use for correcting the data.
#' @param seed A seed is set for reproducibility.
#' @import harmony
#' @details
#' See https://doi.org/10.1038/s41592-019-0619-0 for more details on the Harmony algorithm.
#'
#' @returns fcd with a harmonized expression data frame.
#'
#' @export
harmonize_intensities <- function(fcd,
                                  batch_var,
                                  seed = 91) {

  set.seed(seed)

  harmony_param <- HarmonyMatrix(data_mat = as.matrix(fcd$expr$orig),
                                 meta_data = fcd$anno$cell_anno,
                                 vars_use = batch_var,
                                 do_pca = FALSE)

  fcd[["expr"]][["norm"]] <- as.data.frame(harmony_param)

  return(fcd)

}


#' harmonize_PCA
#'
#' @title harmonize_PCA
#' @description Harmonize the Principal Component Analysis.
#' @param fcd flow cytometry dataset.
#' @param data_slot name of the PCA data slot to use to harmonize. If no prefix was added the, \code(orig).
#' @param batch_var vector of column names from \code{fcd$anno$cell_anno} to use for correcting the data.
#' @param seed A seed is set for reproducibility.
#' @param prefix Prefix for the output.
#' @details
#' See https://doi.org/10.1038/s41592-019-0619-0 for more details on the Harmony algorithm.
#' @return harmonize_PCA
#'
#' @export
harmonize_PCA <- function(fcd,
                          data_slot = "orig",
                          batch_var,
                          seed = 91,
                          prefix = NULL) {

  set.seed(seed)

  if (is.null(prefix)) {

    fcd[["pca"]][["norm"]] <- HarmonyMatrix(data_mat = as.matrix(fcd$pca[[data_slot]]),
                                            meta_data = fcd$anno$cell_anno,
                                            vars_use = batch_var,
                                            do_pca = FALSE)

  } else {

    fcd[["pca"]][[paste(prefix, "norm", sep = "_")]] <- HarmonyMatrix(data_mat = as.matrix(fcd$pca[[data_slot]]),
                                                                      meta_data = fcd$anno$cell_anno,
                                                                      vars_use = batch_var,
                                                                      do_pca = FALSE)

  }

  return(fcd)

}


# 'train_cytonorm
# '
# '@title train_cytonorm
#' @description
#' Wrapper function around 'CytoNorm.train' from the CytoNorm package.
#' @param fcd flow cytometry dataset
#' @param batch_var Column name of batch variable from \code{fcd$anno$cell_anno}.
#' @param remove_param Parameters/markers which should be excluded for learning the batch effect and training the model.
#' @param seed A seed is set for reproducibility.
#' @param files Vector of FCS file names of reference samples which are used for training the model. If files == NULL, all files contained in the fcd are used.
#' @param FCSpath File path to folder where .fcs files contained in the fcd are stored. This parameter does not need to be provided, unless the folder where the .fcs files are stored has changed.
#' @param FlowSOM_param A list of parameters to pass to the FlowSOM algorithm. Default= list(nCells = 5000, xdim = 5, ydim = 5, nClus = 10, scale= FALSE)
#' @returns The function returns a fcd with the trained model saved in extras.
#' @details
#' train_cytonorm' takes a fcd as an input and learns the batch effect of a given batch variable across reference samples provided by the user using the CytoNorm algorithm. This function returns a fcd with the trained model which can be used as input for the \code{\link{run_cytonorm}} function to normalize samples with the trained model.
#' See https://doi.org/10.1002/cyto.a.23904 for more details.
#' @import CytoNorm
#' @export
#'


train_cytonorm <- function(fcd,
                           batch_var,
                           remove_param = NULL,
                           seed,
                           files = NULL,
                           FCSpath = NULL,
                           FlowSOM_param = list(
                             nCells = 5000,
                             xdim = 5,
                             ydim = 5,
                             nClus = 10,
                             scale = FALSE
                           )) {

  # all files are used for training, if no file names are provided
  if(is.null(files)){
    files <-   unique(fcd$anno$cell_anno$expfcs_filename)
  }

  # check FCSpath
  if(is.null(FCSpath)){
    if(!is.null(fcd[["extras"]][["prep_fcd_param"]][["FCSpath"]])){
      #extract fcs_path from extras slot in fcd
      FCSpath <- fcd[["extras"]][["prep_fcd_param"]][["FCSpath"]]
    } else{
      stop(paste0("There is no FCSpath saved in your fcd. Please provide the path to the folder where the FCS files are stored using the 'FCSpath' parameter."))
    }
  }
  #check if FCSpath exists
  if (!dir.exists(FCSpath)) {
    stop(paste0("FCSpath ", FCSpath,  "does not exist"))
  }

  #check if batch variable exists
  if(!batch_var %in% colnames(fcd$anno$cell_anno)){
    stop(paste0(batch_var, " is no column name of fcd$anno$cell_anno"))
  }

  #parameters to remove
  if(!is.null(fcd[["extras"]][["prep_fcd_param"]][["remove_param"]])){
    #extract remove_param from extras slot in fcd
    remove_param_fcd <- fcd[["extras"]][["prep_fcd_param"]][["remove_param"]]
  } else{
    stop("fcd$extras$prep_fcd_param$remove_param is required but is missing in your fcd")
  }

  if(!is.null(remove_param)){
    remove_param <- unique(remove_param, remove_param_fcd)
  }

  # Set seed for reproducibility
  set.seed(seed)

  #select training data
  train_data <-
    fcd$anno$cell_anno[fcd$anno$cell_anno$expfcs_filename %in% files, c("expfcs_filename", batch_var)]
  train_data <-
    train_data %>% distinct(expfcs_filename, .keep_all = TRUE)
  train_data[["path"]] <-
    file.path(FCSpath, train_data$expfcs_filename)

  #get channels
  ff <-flowCore::read.FCS(train_data$path[1], transformation = "linearize", truncate_max_range = FALSE)
  tmp <- ff@parameters@data[, c("name", "desc")]
  tmp <- na.omit(tmp)
  # remove unwanted channels
  if(length(tmp$desc[tmp$desc %in% remove_param])>0){
    keep <- tmp$desc[!tmp$desc %in% remove_param]
    channels <- tmp[tmp$desc %in% keep, "name"]
  } else {
    channels <-  tmp[, "name"]
  }
  transformList <- flowCore::transformList(channels, cytofTransform)
  message("start CytoNorm.train")

  # train model
  model <- CytoNorm::CytoNorm.train(
    files = train_data$path,
    labels = train_data[[batch_var]],
    channels = channels,
    transformList = transformList,
    FlowSOM.params = list(
      nCells = FlowSOM_param[["nCells"]],
      xdim = FlowSOM_param[["xdim"]],
      ydim = FlowSOM_param[["ydim"]],
      nClus = FlowSOM_param[["nClus"]],
      scale = FlowSOM_param[["scale"]]
    ),
    normMethod.train = QuantileNorm.train,
    normParams = list(nQ = 101,
                      goal = "mean"),
    seed = seed,
    recompute = TRUE,
    truncate_max_range= FALSE
  )
  #add model to fcd
  fcd[["extras"]][["cytonorm_model"]] <- model

  return(fcd)
}


#' run_cytonorm
#'
#' @title run_cytonorm
#' @description
#' Wrapper function around CytoNorm.normalize from the CytoNorm package.
#' @param fcd flow cytometry dataset
#' @param batch_var Column name of batch variable from \code{fcd$anno$cell_anno}.
#' @param keep_fcs Boolean whether to keep the normalized FCS files in \code{output_dir}.
#' @param output_dir Directory to save normalized FCS files temporary or permanently, if keep_fcs == TRUE.
#' @param files Vector of fcs file names of samples which should be normalized. By default all files contained in the flow cytometry dataset are used.
#' @param FCSpath File path to folder where .fcs files contained in the fcd are stored.
#' @param anno_table Path to the annotation table file.
#' @returns fcd with a normalized expression data frame.
#' @details
#' This function assumes that your fcd contains a trained model computed by \code{\link{train.cytonorm}}. The function performs normalization of the samples contained in your fcd. The normalized expression values are added to your fcd and by default FCS files with the normalized values are written to the  \code{output_dir}.
#' @import CytoNorm
#' @export
run_cytonorm <- function(fcd,
                         batch_var,
                         keep_fcs = TRUE,
                         output_dir = paste0("./CytoNorm_output_", Sys.Date()),
                         files= NULL,
                         FCSpath=NULL,
                         anno_table= NULL){

  # all files are used for normalization, if no files names are provided
  if(is.null(files)){
    files <-   unique(fcd$anno$cell_anno$expfcs_filename)
  }

  # check FCSpath
  if(is.null(FCSpath)){
    if(!is.null(fcd[["extras"]][["prep_fcd_param"]][["FCSpath"]])){
      #extract fcs_path from extras slot in fcd
      FCSpath <- fcd[["extras"]][["prep_fcd_param"]][["FCSpath"]]
    } else{
      stop(paste0("There is no FCSpath saved in your fcd. Please provide the path to the folder where the FCS files are stored using the 'FCSpath' parameter."))
    }
  }
  #check if FCSpath exists
  if (!dir.exists(FCSpath)) {
    stop("'FCSpath' does not exist.")
  }
  #check if fcd contains model
  if (!"cytonorm_model" %in% names(fcd[["extras"]])) {
    stop("Your fcd does not contain the cytonorm_model.")
  }

  #select data to normalize
  data <-
    fcd$anno$cell_anno[fcd$anno$cell_anno$expfcs_filename %in% files, c("expfcs_filename", batch_var)]
  data <- data %>% distinct(expfcs_filename, .keep_all = TRUE)
  data[["path"]] <-file.path(FCSpath, data$expfcs_filename)

  # extract channels used for calculating the model
  channels <- fcd[["extras"]][["cytonorm_model"]][["clusterRes"]][["1"]][["channels"]]
  transformList <- flowCore::transformList(channels, CytoNorm::cytofTransform)
  transformList.reverse <- flowCore::transformList(channels, CytoNorm::cytofTransform.reverse)

  message("start normalization")
  # normalize fcs files
  CytoNorm::CytoNorm.normalize(model = fcd[["extras"]][["cytonorm_model"]],
                               files = data$path,
                               labels = data[[batch_var]],
                               transformList = transformList,
                               transformList.reverse = transformList.reverse,
                               normMethod.normalize = QuantileNorm.normalize,
                               outputDir = output_dir,
                               prefix= "", # no prefix to use same fcd annotation table
                               clean = TRUE,
                               verbose = FALSE,
                               truncate_max_range= FALSE)

  #prep_fcd_params
  if(!is.null(fcd[["extras"]][["prep_fcd_param"]])){
    #check if all required parameters are saved in fcd
    param_names <- c("ceil", "transformation", "remove_param", "filename_col", "seed", "separator_anno")
    if(sum(param_names %in% names(fcd[["extras"]][["prep_fcd_param"]])) != length(param_names
    )){
      stop(paste0(setdiff(param_names, names(fcd[["extras"]][["prep_fcd_param"]])), " is missing in fcd$extras$prep_fcd_param."))
    }
    #extract prep_fcd_param from extras slot in fcd
    prep_fcd_param <- list(ceil = fcd[["extras"]][["prep_fcd_param"]][["ceil"]],
                           transformation =fcd[["extras"]][["prep_fcd_param"]][["transformation"]],
                           remove_param= fcd[["extras"]][["prep_fcd_param"]][["remove_param"]],
                           # anno_table = fcd[["extras"]][["prep_fcd_param"]][["anno_table"]],
                           filename_col= fcd[["extras"]][["prep_fcd_param"]][["filename_col"]],
                           seed= fcd[["extras"]][["prep_fcd_param"]][["seed"]],
                           separator_anno = fcd[["extras"]][["prep_fcd_param"]][["separator_anno"]])
    #check if anno_table parameter was provided
    if(!is.null(anno_table)){
      prep_fcd_param[["anno_table"]] <- anno_table
    } else{
      if(!is.null(fcd[["extras"]][["prep_fcd_param"]][["anno_table"]])){
        prep_fcd_param[["anno_table"]] <-fcd[["extras"]][["prep_fcd_param"]][["anno_table"]]
      }else{
        stop(paste0("There is no anno_table parameter saved in your fcd. Please provide the path to the folder where the annotation file using the 'anno_table' parameter."))
      }
    }
    #check if anno table file exists
    if (!file.exists( prep_fcd_param[["anno_table"]])) {
      stop(paste0("The annotation file ", prep_fcd_param[["anno_table"]], " cannot be found. Please provide the correct path to the annotation file using the 'anno_table' parameter."))
    }
  } else{
    stop("The fcd$extras$prep_fcd_param is required but cannot be found  your fcd.")
  }


  # read in normalized FCS files
  message("adding normalized expression data to fcd")
  fcd_norm<- prep_fcd(FCSpath = output_dir,
                      ceil = prep_fcd_param[["ceil"]],
                      useCSV = FALSE,
                      transformation = prep_fcd_param[["transformation"]],
                      remove_param = prep_fcd_param[["remove_param"]],
                      anno_table = prep_fcd_param[["anno_table"]],
                      filename_col = prep_fcd_param[["filename_col"]],
                      seed = prep_fcd_param[["seed"]])


  #add normalized expression data frame to original fcd
  if(identical(colnames(fcd$expr$orig), colnames(fcd_norm$expr$orig)))
  {if(identical(rownames(fcd$expr$orig), rownames(fcd_norm$expr$orig)))
    #add normalized expression data frame
  { fcd[["expr"]][["norm"]] <- fcd_norm[["expr"]][["orig"]]
  }else{
    stop("The rownames of the original and the normalized expression data frame are not identical.")
  }
  } else{
    stop("The normalized expression table cannot be added to the fcd because the colnames of the original and normalized expression data frame are not identical.")
  }


  #cleaning
  if(keep_fcs== FALSE){
    message("removing temporary fcs files")
    tmp_files <- file.path(output_dir,list.files(output_dir))
    file.remove(tmp_files)
    file.remove(list.dirs(output_dir), recursive= TRUE)
  }

  return(fcd)
}
