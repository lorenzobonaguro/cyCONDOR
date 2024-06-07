#' harmonize_intensities
#'
#' @title harmonize_intensities
#' @description Harmonize the expression values.
#' @param fcd flow cytometry dataset.
#' @param batch vector of column names to use for correcting the data.
#' @param seed Seed used for the randomization steps.
#' @import harmony
#' @return harmonize_intensities
#'
#' @export
harmonize_intensities <- function(fcd, batch, seed) {

  set.seed(seed)

  harmony_param <- HarmonyMatrix(data_mat = as.matrix(fcd$expr$orig),
                                 meta_data = fcd$anno$cell_anno,
                                 vars_use = batch,
                                 do_pca = FALSE)

  fcd[["expr"]][["norm"]] <- as.data.frame(harmony_param)

  return(fcd)

}


#' harmonize_PCA
#'
#' @title harmonize_PCA
#' @description Harmonize the Principal Component Analysis.
#' @param fcd flow cytometry dataset.
#' @param data_slot name of the PCA data slot to use to harmonize. If no prefix was added the, *orig*.
#' @param batch vector of column names to use for correcting the data.
#' @param seed Seed used for the randomization steps.
#' @param prefix Prefix for the output.
#' @return harmonize_PCA
#'
#' @export
harmonize_PCA <- function(fcd, data_slot = "orig", batch, seed, prefix = NULL) {

  set.seed(seed)

  if (is.null(prefix)) {

    fcd[["pca"]][["norm"]] <- HarmonyMatrix(data_mat = as.matrix(fcd$pca[[data_slot]]),
                                            meta_data = fcd$anno$cell_anno,
                                            vars_use = batch,
                                            do_pca = FALSE)

  } else {

    fcd[["pca"]][[paste(prefix, "norm", sep = "_")]] <- HarmonyMatrix(data_mat = as.matrix(fcd$pca[[data_slot]]),
                                                                      meta_data = fcd$anno$cell_anno,
                                                                      vars_use = batch,
                                                                      do_pca = FALSE)

  }

  return(fcd)

}


#'train_cytonorm
#'
#'@title train_cytonorm
#'@description
#'wrapper function around CytoNorm.train from the CytoNorm package.
#'@param fcd flow cytometry dataset
#'@param files Vector of fcs file names of reference samples which are used for training the model. If files == NULL, all files contained in the flow cytometry dataset are used.
#'@param FCSpath File path to folder where .fcs files contained in the fcd are stored.
#'@param batch_var Column name of batch variable from fcd$anno$cell_anno
#'@param remove_param Parameters which should be excluded for normalization
#'@param FlowSOM_param nClus is important to set, refer to help function
#'@returns The function returns a fcd with the trained model saved as extras.
#'@details
#'train_cytonorm' takes a fcd as an input and learns the batch effect of a given batch variable across reference samples provided by the user using CytoNorm. This function returns a fcd with the trained model which can then be used as input for the 'run_cytonorm' function to normalize samples with the trained model.
#'See https://doi.org/10.1002/cyto.a.23904 for more details.
#'@import CytoNorm
#'@import flowCore


train_cytonorm <- function(fcd,
                           files = NULL,
                           FCSpath,
                           batch_var,
                           remove_param = NULL,
                           FlowSOM_param = list(
                             nCells = 5000,
                             xdim = 5,
                             ydim = 5,
                             nClus = 10,
                             scale = FALSE
                           ),
                           seed) {

  # all files are used for training, if no files names are provided
  if(is.null(files)){
    files <-   unique(fcd$anno$cell_anno$expfcs_filename)
  }

  #check if FCS path exists
  if (!dir.exists(FCSpath)) {
    stop("FCSpath does not exist")
  }
  #check if batch variable exists
  if(!batch_var %in% colnames(fcd$anno$cell_anno)){
    stop(paste0(batch_var, " is no column name of fcd$anno$cell_anno"))
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


#'run_cytonorm
#'
#'@title run_cytonorm
#'@description
#'@param fcd flow cytometry dataset
#'@param files Vector of fcs file names of samples which should be normalized. By default all files contained in the flow cytometry dataset are used.
#'@param FCSpath File path to folder where .fcs files contained in the fcd are stored.
#'@param batch_var Column name of batch variable from fcd$anno$cell_anno
#'@param output_dir Directory to save normalized fcs files.
#'@param prep_fcd_param Parameters for prep_fcd function. The same parameters as for the unnormalized fcd have to be provided.
#'@param keep_fcs Boolean whether to keep the normalized FCS files in output_dir.
#'@returns fcd with a normalized expression data frame.
#'@details
#'to be added
#'@import CytoNorm
#'@import flowcore
run_cytonorm <- function(fcd,
                         files= NULL,
                         FCSpath,
                         batch_var,
                         output_dir = paste0("./CytoNorm_output_", Sys.Date()),
                         prep_fcd_param= list(ceil,
                                              transformation,
                                              remove_param,
                                              anno_table,
                                              filename_col,
                                              seed),
                         keep_fcs = TRUE){

  # all files are used for normalization, if no files names are provided
  if(is.null(files)){
    files <-   unique(fcd$anno$cell_anno$expfcs_filename)
  }

  #check if FCS path exists
  if (!dir.exists(FCSpath)) {
    stop("FCSpath does not exist")
  }
  #check if fcd contains model
  if (!"cytonorm_model" %in% names(fcd[["extras"]])) {
    stop("fcd does not contain cytonorm_model")
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
  CytoNorm::CytoNorm.normalize(model = fcd[["extras"]][["cytonorm_model"]], # suppress warning
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

  # read in normalized FCS files
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
    #add cytonorm expression data frame
  { fcd[["expr"]][["norm"]] <- fcd_norm[["expr"]][["orig"]]
  }else{
    stop("The rownames of the original and the normalized expression data frame are not identical")
  }
  } else{
    stop("The normalized expression table cannot be added to fcd because the colnames of the original and normalized expression data frame are not identical")
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
