#' crt_transform
#'
#' @title crt_transform
#' @description Data transformation, this function run within the prep_fcd wrapper.
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
#' @description Data transformation, this function run within the prep_fcd wrapper, the logicle tranformation are derived from Cytofkit.
#' @param keep Vector of the parameter to keep in the analysis.
#' @param original_data Original data
#' @param transformation transformation to perform.
#' @return transformed flow cytometry dataset
#'
#' @export
transform_data <- function(keep, transformation, original_data){

  # Save temp df
  transf_data <- original_data

  # Define fixed transformations
  CyTOFlgcl <- logicleTransform(w=0.25, t=16409, m=4.5, a=0)

  # Transform the data
  for(paramName in as.character(keep)){

    if(transformation == "cytof_logi"){
      dataNum <- which(colnames(original_data)==paramName)
      temp <- apply(original_data[,dataNum,drop=F],2, CyTOFlgcl)
      transf_data[,dataNum] <- temp
    }

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
#' @description Wrapping function to prepare a flow cytometry dataset
#' @param data_path Folder where the .fcs files or .csv files are stored.
#' @param max_cell Number of cells to use for each file (set to a high number if you want to use all available events).
#' @param useCSV Flag if the input are .csv files and not .fcs (experimental).
#' @param transformation Transformation to perform.
#' @param remove_param Parameters to remove from the transformation, "inTime" should be kept.
#' @param anno_table path to the annotation table file.
#' @param filename_col Name of the column containing the file name matching with the .fcs files.
#' @param seed A seed is set for reproducibility.
#' @param separator_anno separator used in the annotation file.
#' @param separator_fc_csv separator used in the fc csv files.
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
  if (!transformation %in% c("cytof_logi", "clr", "arcsinh", "auto_logi")) {
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
#' @description read_flowjo_workspace and prepare the condor object
#' @param data_gs Gate Set object from flowWorkspace Package.
#' @param pop Gate to keep for downstream analysis (default: 'root').
#' @param gate_list Gate List of the FlowJo Workspace.
#' @param inverse.transform Logical: if the data should be reverse transformed of kept with FlowJo transformation (default = FALSE).
#' @param transformation If inverse.transform = TRUE, type of new transformation to perform (see nfTransform).
#' @param remove_param Parameters to be removed from the condor object.
#' @param merge_anno Logical: If sample anno should be merged to the condor object.
#' @param anno_table Path to annotation table.
#' @param separator_anno Separator of the .csv annotation table.
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
                     anno_table,
                     separator_anno) {

  fs <- flowWorkspace::gs_pop_get_data(obj = data_gs, y = "root", inverse.transform = inverse.transform) %>% flowWorkspace::cytoset_to_flowSet()

  filenames <- rownames(fs@phenoData)

  NumBC <- length(fs)

  FFdata <- NULL

  OrigNames <-fs[[1]]@parameters$name

  for (FFs in 1:NumBC){
    FFa <- exprs(fs[[FFs]])

    #Fixup column names
    colnames(FFa) <- fs[[FFs]]@parameters$desc
    empties <- which(is.na(colnames(FFa)) | colnames(FFa)== " ")
    colnames(FFa)[empties] <- fs[[FFs]]@parameters$name[empties]
    fs[[FFs]]@parameters$desc <- colnames(FFa)

    #Add file label
    FFa <- cbind(FFa,rep(FFs,dim(FFa)[1]))
    colnames(FFa)[dim(FFa)[2]] <- "InFile"

    # Add the gating info
    for (gate in gate_list) {

      FFa <- cbind(FFa, gh_pop_get_indices(gs[[FFs]], y = gate))
      colnames(FFa)[dim(FFa)[2]] <- gate

    }

    #Concatenate
    FFdata <- rbind(FFdata,FFa)
  }

  FFdata <- as.data.frame(FFdata)

  if (inverse.transform == TRUE) {

    ## Data Transformation
    keep <- colnames(raw_data)[!colnames(raw_data) %in% remove_param]

    raw_data <- raw_data[,which(colnames(raw_data) %in% keep)]

    ## Check if transformation is a valid value
    if (!transformation %in% c("cytof_logi", "fluor_logi", "clr", "arcsinh", "auto_logi")) {
     stop(paste0(transformation, " is not a valid transformation method"))
    }

    trans_data <- transform_data(keep = keep, transformation = transformation, original_data = raw_data)

    ## Clean the dataframe
    df <- cbind(trans_data, expfcs_filename=data$merged_df[,"InFile"])
    df <- as.data.frame(df)
    df$expfcs_filename <- as.factor(df$expfcs_filename)
    df$expfcs_filename <- factor(df$expfcs_filename, labels = data$FcsFiles)

  } else {

    df <- FFdata
    colnames(df)[colnames(df) == "InFile"] <- "expfcs_filename"
    df$expfcs_filename <- as.factor(df$expfcs_filename)
    df$expfcs_filename <- factor(df$expfcs_filename, labels = data$FcsFiles)

  }

  if (merge_anno == TRUE) {c

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

  class(fcd) <- "flow_cytometry_dataframe"

  return(fcd)

}

