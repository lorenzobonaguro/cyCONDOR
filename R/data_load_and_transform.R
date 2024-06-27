#' crt_transform
#'
#' @title crt_transform
#' @description Data transformation, this function runs within the \code{\link{prep_fcd}} wrapper.
#' @param x matrix to transform
#' @return crt_transform
#'
#' @export
clr <- function(x) {
  return(log1p(x = x/(exp(x = sum(log1p(x = x[x > 0]), na.rm = TRUE)/length(x = x)))))
  }

#' read_data
#'
#' @title read_data
#' @description Load .fcs or .csv files into a dataframe and prepare the condor object.
#' @param FCSpath Path to the .fcs files.
#' @param max_cells number of cells to subset.
#' @param useCSV Logical, if input is .csv and not .fcs.
#' @param separator Separator used the flow csv files (if loading from csv).
#' @import flowCore
#' @import reshape2
#' @import dplyr
#' @import cowplot
#' @return load flow cytometry dataset
#'
#' @export
read_data <- function(FCSpath, max_cells, useCSV, separator){

  if(useCSV == FALSE){

    # Read FCS files

    FcsFiles <- list.files(path = FCSpath, pattern = ".fcs")
    flow_set = list()
    for(FileNum in 1:length(FcsFiles)){
      flow_set[[FileNum]] <- read.FCS(paste0(FCSpath,"/",FcsFiles[FileNum]),
                                      transformation =F,
                                      ignore.text.offset=T,
                                      truncate_max_range=F,
                                      emptyValue = F)}

    # Arrange the dataset
    num_files <- length(flow_set)
    merged_df <- NULL

    for (single_file in 1:num_files){
      single_file_red <- exprs(flow_set[[single_file]])

      ## Downsample if needed
      if (nrow(single_file_red)<=max_cells)
        single_file_red <- single_file_red
      else
        single_file_red <- single_file_red[sample(nrow(single_file_red),max_cells,replace=F),]

      #Adjust colnames
      colnames(single_file_red) <- flow_set[[single_file]]@parameters$desc
      colnames(single_file_red)[which(is.na(colnames(single_file_red)) | colnames(single_file_red)== " ")] <- flow_set[[single_file]]@parameters$name[which(is.na(colnames(single_file_red)) | colnames(single_file_red)== " ")]

      #Include column with file label for tracking
      single_file_red <- cbind(single_file_red,rep(single_file,dim(single_file_red)[1]))
      colnames(single_file_red)[dim(single_file_red)[2]] <- "InFile"

      #Merge
      merged_df <- rbind(merged_df,single_file_red)

    }} else {

      csvfilenames <- list.files(path = FCSpath, pattern="*.csv")
      FcsFiles <- csvfilenames
      csvdata <- lapply(paste0(FCSpath,"//",csvfilenames),function(x) read.delim(x, check.names = F, sep = separator))
      num_files <- length(csvdata)
      merged_df<-NULL
      for (single_file in 1:num_files){
        single_file_red <- csvdata[[single_file]]

        ## Downsample if needed
        if (nrow(single_file_red)<=max_cells)
          single_file_red <- single_file_red
        else
          single_file_red <- single_file_red[sample(nrow(single_file_red),max_cells,replace=F),]

        single_file_red <- cbind(single_file_red,rep(single_file,dim(single_file_red)[1]))
        colnames(single_file_red)[dim(single_file_red)[2]] <- "InFile"

        #Merge
        merged_df <- rbind(merged_df,single_file_red)
      }
    }

  return(list(merged_df=merged_df, FcsFiles = FcsFiles))
}

#' transform_data
#'
#' @title transform_data
#' @description Data transformation, this function runs within the \code{\link{prep_fcd}} wrapper, the logicle tranformation are derived from Cytofkit.
#' @param keep Vector of the parameter to keep in the analysis.
#' @param original_data Original data
#' @param transformation transformation to perform.
#' @return transformed flow cytometry dataset
#'
#' @export
transform_data <- function(keep, transformation, original_data){

  # Save temp df
  transf_data <- original_data

  # Transform the data
  for(paramName in as.character(keep)){

    if(transformation == "clr" ){
      dataNum <- which(colnames(original_data)==paramName)
      temp <- apply(original_data[,dataNum,drop=F],2, clr)
      transf_data[,dataNum] <- temp
    }

    if(transformation == "arcsinh" ){
      dataNum <- which(colnames(original_data)==paramName)
      temp <- original_data[,dataNum,drop=F] / 5
      temp <- asinh(temp)
      transf_data[,dataNum] <- temp
    }

    if(transformation == "auto_logi"){
      q<-0.05
      m<-4.5
      d <- original_data[,paramName]
      w <- 0
      t <- max(d)
      nd <- d[d < 0]
      nThres <- quantile(nd, 0.25) - 1.5 * IQR(nd)
      nd <- nd[nd >= nThres]
      if (length(nd)) {
        r <- .Machine$double.eps + quantile(nd, q)
        if (10^m * abs(r) <= t) {
          w <- 0
        }
        else {
          w <- (m - log10(t/abs(r)))/2
          if (is.nan(w) || w > 2) {
            warning(paste0("autoLgcl failed for channel: ",
                           paramName, "; using default fluor logicle transformation!"))
            w <- 0.1
            t <- 500000
            m <- 4.5
          }
        }
      }
      templgcl <- logicleTransform(w=w, t=t, m=4.5, a=0)
      dataNum <- which(colnames(original_data)==paramName)
      temp <- apply(original_data[,dataNum,drop=F],2, templgcl)
      transf_data[,dataNum] <- temp
      print(paste0(paramName, " w= ",w," t= ",t))
    }

  }
  return(transf_data)
}


#' prep_fcd
#'
#' @title prep_fcd
#' @description Loading and transforming the data to create a flow cytometry dataset from FCS files for the analysis with the cyCONDOR workflow.
#' @param data_path Folder where the .fcs files or .csv files are stored.
#' @param max_cell Number of cells to use for each file (set to a high number if you want to use all available events).
#' @param useCSV Flag if the input are .csv files and not .fcs (experimental).
#' @param transformation Transformation to perform. Select one of the following: \code{"auto_log"} (autologicle, recommended for flow cytometry data), \code{"arcsinh"} (arcsinh transformation with cofactor 5), \code{"clr"} (centered-log-ratio) or \code{"none"} (no transformation).
#' @param remove_param Parameters to remove from the transformation, "inTime" should be kept. The selected parameters will not be included in the \code{fcd}.
#' @param anno_table Path to the annotation table text file. The annotation table should contain one column with the file names of all .fcs or .csv files to read in and optionally additional columns with further sample information (e.g. "sample_id", "condition").
#' @param filename_col Name of the column of the \code{anno_table} containing the file name matching with the .fcs/.csv files.
#' @param seed A seed is set for reproducibility.
#' @param separator_anno Separator used in the annotation file, by default \code{separator_anno = ","}.
#' @param separator_fc_csv Separator used in the cytometry data .csv files, by default \code{separator_anno = ","}.
#' @details The \code{prep_fcd} is a wrapper function to read in the files, subset to \code{max_cell}, transform the data and create a 'flow cytometry dataframe' (\code{fcd}).
#' @return An object of class 'flow cytometry dataframe' (\code{fcd}) is returned.
#' @import readr
#' @import readxl
#' @import stringr
#' @import Rmisc
#' @return prep_fcd
#'
#' @export
prep_fcd <- function(data_path,
                     max_cell,
                     useCSV = FALSE,
                     transformation,
                     remove_param,
                     anno_table,
                     filename_col,
                     seed = 91,
                     separator_anno = ",",
                     separator_fc_csv = ",") {

  # Set seed for reproducibility
  set.seed(seed)

  ## Load the data
  data <- read_data(FCSpath = data_path, max_cells = max_cell, useCSV = useCSV, separator = separator_fc_csv)

  raw_data <- as.matrix(data$merged_df) # Take the dataframe with the intensity values

  ## Data Transformation
  keep <- colnames(raw_data)[!colnames(raw_data) %in% remove_param]

  raw_data <- raw_data[,which(colnames(raw_data) %in% keep)]

  ## Check if transformation is a valid value
  if (!transformation %in% c("clr", "arcsinh", "auto_logi", "none")) {
    stop(paste0(transformation, " is not a valid transformation method"))
  }

  trans_data <- transform_data(keep = keep, transformation = transformation, original_data = raw_data)

  ## Clean the dataframe
  df <- cbind(trans_data, expfcs_filename=data$merged_df[,"InFile"])
  df <- as.data.frame(df)
  df$expfcs_filename <- as.factor(df$expfcs_filename)
  df$expfcs_filename <- factor(df$expfcs_filename, labels = data$FcsFiles)

  ## Now add the annotation (as csv file)
  anno <- read.delim(anno_table, sep = separator_anno)

  df <- merge(df, anno, by.x = "expfcs_filename", by.y = filename_col)

  ## Give the unique rownames
  rownames(df) <- paste(df$expfcs_filename, rownames(df), sep = "_")

  ## Prepare the final object
  fcd <- list()

  fcd[["expr"]][["orig"]] <- df[ ,colnames(df) %in% keep]
  fcd[["anno"]][["cell_anno"]] <- df[ ,!colnames(df) %in% keep]

  #save import parameters
  fcd[["extras"]][["prep_fcd_param"]] <- list(data_path = data_path,
                                              max_cell = max_cell,
                                              transformation = transformation,
                                              remove_param= remove_param,
                                              anno_table = anno_table,
                                              filename_col= filename_col,
                                              seed= seed,
                                              separator_anno = separator_anno,
                                              separator_fc_csv = separator_fc_csv)

  class(fcd) <- "flow_cytometry_dataframe"

  return(fcd)

}


#' Read FlowJo workspace
#'
#' @title Read FlowJo Workspace
#' @description Loading and transforming the data to create a flow cytometry dataset from a Gate Set object for the analysis with the cyCONDOR workflow.
#' @param data_gs Gate Set object, e.g. created by using \code{\link[CytoML]{open_flowjo_xml}} and \code{\link[CytoML]{flowjo_to_gatingset}} from the CytoML package.
#' @param pop Gate to keep for downstream analysis (default: 'root').
#' @param gate_list Gate List of the FlowJo Workspace.
#' @param inverse.transform Logical: if the data should be reverse transformed of kept with FlowJo transformation (default = FALSE).
#' @param transformation If \code{inverse.transform = TRUE}, type of new transformation to perform. Select one of the following: \code{"auto_log"} (autologicle, recommended for flow cytometry data), \code{"arcsinh"} (arcsinh transformation with cofactor 5), \code{"clr"} (centered-log-ratio) or \code{"none"} (no transformation).
#' @param remove_param Parameters to be removed from the condor object.
#' @param merge_anno Logical: If sample anno should be merged to the \code{fcd}.
#' @param anno_table If \code{merge_anno = TRUE}, path to the annotation table text file. The annotation table should contain one column with the file names of all .fcs or .csv files to read in and optionally additional columns with further sample information (e.g. "sample_id", "condition").
#' @param separator_anno Separator used in the annotation file, by default \code{separator_anno = ","}.
#' @param filename_col Name of the column of the \code{anno_table} containing the file name matching with the .fcs files.
#' @return An object of class 'flow cytometry dataframe' (\code{fcd}) is returned.
#' @import flowWorkspace
#' @import Biobase
#' @import CytoML
#' @return read_flowjo_workspace
#'
#' @export
prep_fjw <- function(data_gs,
                     pop = "root",
                     gate_list = nodelist,
                     inverse.transform = FALSE,
                     transformation,
                     remove_param,
                     merge_anno = FALSE,
                     anno_table = NULL,
                     separator_anno = ",",
                     filename_col = NULL) {

  fs <- flowWorkspace::gs_pop_get_data(obj = data_gs, y = "root", inverse.transform = inverse.transform) %>% flowWorkspace::cytoset_to_flowSet()

  filenames <- rownames(fs@phenoData)

  num_files <- length(fs)

  data <- NULL

  for (single_file in 1:num_files){
    raw_file <- exprs(fs[[single_file]])

    #Fixup column names
    colnames(raw_file) <- fs[[single_file]]@parameters$desc
    colnames(raw_file)[which(is.na(colnames(raw_file)) | colnames(raw_file)== " ")] <- fs[[single_file]]@parameters$name[which(is.na(colnames(raw_file)) | colnames(raw_file)== " ")]
    fs[[single_file]]@parameters$desc <- colnames(raw_file)

    #Add file label
    raw_file <- cbind(raw_file,rep(single_file,dim(raw_file)[1]))
    colnames(raw_file)[dim(raw_file)[2]] <- "InFile"

    # Add the gating info
    for (gate in gate_list) {

      raw_file <- cbind(raw_file, gh_pop_get_indices(gs[[single_file]], y = gate))
      colnames(raw_file)[dim(raw_file)[2]] <- gate

    }

    #Concatenate
    data <- rbind(data,raw_file)
  }

  raw_data <- as.data.frame(data)

  if (inverse.transform == TRUE) {

    ## Data Transformation
    keep <- fs[[1]]@parameters$desc[!fs[[1]]@parameters$desc %in% remove_param]

    raw_data <- raw_data[,which(colnames(raw_data) %in% keep)]

    ## Check if transformation is a valid value
    if (!transformation %in% c("clr", "arcsinh", "auto_logi", "none")) {
     stop(paste0(transformation, " is not a valid transformation method"))
    }

    trans_data <- transform_data(keep = keep, transformation = transformation, original_data = raw_data)

    ## Clean the dataframe
    df <- cbind(trans_data, data[, !colnames(data) %in% fs[[1]]@parameters$desc])
    df <- as.data.frame(df)
    colnames(df)[colnames(df) == "InFile"] <- "expfcs_filename"
    df$expfcs_filename <- as.factor(df$expfcs_filename)
    df$expfcs_filename <- factor(df$expfcs_filename, labels = filenames)

  } else {

    df <- data
    colnames(df)[colnames(df) == "InFile"] <- "expfcs_filename"
    df$expfcs_filename <- as.factor(df$expfcs_filename)
    df$expfcs_filename <- factor(df$expfcs_filename, labels = filenames)

  }

  if (merge_anno == TRUE) {

    ## Now add the annotation (as csv file)
    anno <- read.delim(anno_table, sep = separator_anno)

    df <- merge(df, anno, by.x = "expfcs_filename", by.y = filename_col)

  }

  ## Give the unique rownames
  rownames(df) <- paste(df$expfcs_filename, rownames(df), sep = "_")

  ## Prepare the final object
  fcd <- list()

  fcd[["expr"]][["orig"]] <- df[ ,colnames(df) %in% fs[[1]]@parameters$desc]
  fcd[["anno"]][["cell_anno"]] <- df[ ,!colnames(df) %in% fs[[1]]@parameters$desc]

  #save import parameters
  fcd[["extras"]][["prep_fcd_param"]] <- list(#data_path = data_path,
                                             #max_cell = max_cell,
                                              transformation = transformation,
                                              remove_param= remove_param,
                                              anno_table = anno_table,
                                              filename_col= filename_col,
                                             #seed= seed,
                                              separator_anno = separator_anno)


  class(fcd) <- "flow_cytometry_dataframe"

  return(fcd)

}

