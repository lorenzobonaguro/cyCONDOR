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
#' @param data_path Path to the .fcs or .csv files.
#' @param max_cells number of cells to subset.
#' @param useCSV Logical, if input is .csv and not .fcs.
#' @param separator Separator used the flow csv files (if loading from csv).
#' @param simple_names If TRUE only the channel description is used to name the column, if FALSE both channel name and description are pasted together.
#' @param truncate_max_range From FlowCore: logical type. Default is FALSE. can be optionally turned off to avoid truncating the extreme positive value to the instrument measurement range .i.e.'$PnR'.
#' @param emptyValue From FlowCore: boolean indicating whether or not we allow empty value for keyword values in TEXT segment. It affects how the double delimiters are treated. IF TRUE, The double delimiters are parsed as a pair of start and end single delimiter for an empty value. Otherwise, double delimiters are parsed one part of string as the keyword value. default is TRUE.
#' @param ignore.text.offset From FlowCore: whether to ignore the keyword values in TEXT segment when they don't agree with the HEADER. Default is FALSE, which throws the error when such discrepancy is found. User can turn it on to ignore TEXT segment when he is sure of the accuracy of HEADER so that the file still can be read.
#' @param verbose Default FALSE, if TRUE the at each file loaded something is printed in the screen.
#' @param cross_path_with_anno Defautl FALSE. If TRUE is the 'data_path' contains more files then the annotation table only the overlap will be loaded.
#' @param anno_table Passed from 'prep_fcd'
#' @param separator_anno Passed from 'prep_fcd'
#' @param filename_col Passed from 'prep_fcs'
#' @import flowCore
#' @import reshape2
#' @import dplyr
#' @import cowplot
#' @return load flow cytometry dataset
#'
#' @export
read_data <- function(data_path,
                      max_cells,
                      useCSV,
                      separator,
                      simple_names,
                      truncate_max_range,
                      emptyValue,
                      ignore.text.offset,
                      verbose,
                      cross_path_with_anno,
                      anno_table,
                      separator_anno,
                      filename_col){

  if(useCSV == FALSE){

    # Read FCS files

    data_files <- list.files(path = data_path, pattern = ".fcs")

    if (cross_path_with_anno == TRUE) {

      anno_files <- read.delim(anno_table, sep = separator_anno)[[filename_col]]

      data_files <- data_files[data_files %in% anno_files]

    }

    merged_df <- NULL
    for(FileNum in 1:length(data_files)){

      if (verbose == TRUE) {

        print(paste0("Loading file ", FileNum, " out of ", length(data_files)))

      }

      flow_frame_single <- read.FCS(paste0(data_path,"/",data_files[FileNum]),
                                    transformation =F,
                                    ignore.text.offset=ignore.text.offset,
                                    truncate_max_range=truncate_max_range,
                                    emptyValue = emptyValue)

      single_file_red <- exprs(flow_frame_single)

      if (dim(single_file_red)[1] == 0) {

        stop(paste0("File: ", data_files[FileNum], " contains no valid event"))

      }

      ## Downsample if needed
      if (nrow(single_file_red) <= max_cells) {
        single_file_red <- single_file_red
      } else {
        single_file_red <- single_file_red[sample(nrow(single_file_red),max_cells,replace=F),]
      }

      #Adjust colnames
      if (simple_names == TRUE) {

        colnames(single_file_red) <- flow_frame_single@parameters@data$desc
        colnames(single_file_red)[which(is.na(colnames(single_file_red)) | colnames(single_file_red)== " ")] <- flow_frame_single@parameters@data$name[which(is.na(colnames(single_file_red)) | colnames(single_file_red)== " ")]

      } else {

        colnames(single_file_red) <- ifelse(is.na(flow_frame_single@parameters@data$desc),
                                            flow_frame_single@parameters@data$name,
                                            paste0(flow_frame_single@parameters@data$name, "_", flow_frame_single@parameters@data$desc))

      }


      #Include column with file index for tracking and later merging annotation
      single_file_red <- cbind(single_file_red,rep(FileNum,dim(single_file_red)[1]))
      colnames(single_file_red)[dim(single_file_red)[2]] <- "InFile"

      #Merge into a single data.frame
      merged_df <- rbind(merged_df,single_file_red)

    }
  } else { # read csv

    data_files <- list.files(path = data_path, pattern = ".csv")

    if (cross_path_with_anno == TRUE) {

      anno_files <- read.delim(anno_table, sep = separator_anno)[[filename_col]]

      data_files <- data_files[data_files %in% anno_files]

    }

    merged_df <- NULL

    for (FileNum in 1:length(data_files)){

      if (verbose == TRUE) {

        print(paste0("Loading file ", FileNum, " out of ", length(data_files)))

      }

      single_file_red <- read.delim(paste0(data_path,"/",data_files[FileNum]), check.names = F, sep = separator)

      if (dim(single_file_red)[1] == 0) {

        stop(paste0("File: ", data_files[FileNum], " contains no valid event"))

      }

      ## Downsample if needed
      if (nrow(single_file_red) <= max_cells){
        single_file_red <- single_file_red
      } else {
        single_file_red <- single_file_red[sample(nrow(single_file_red),max_cells,replace=F),]
      }

      #Include column with file index for tracking and later merging annotation
      single_file_red <- cbind(single_file_red,rep(FileNum,dim(single_file_red)[1]))
      colnames(single_file_red)[dim(single_file_red)[2]] <- "InFile"

      #Merge
      merged_df <- rbind(merged_df,single_file_red)
    }
  }

  return(list(merged_df=merged_df, data_files = data_files))
}

#' transform_data
#'
#' @title transform_data
#' @description Data transformation, this function runs within the \code{\link{prep_fcd}} wrapper, the logicle tranformation are derived from Cytofkit.
#' @param keep Vector of the parameter to keep in the analysis.
#' @param original_data Original data
#' @param transformation transformation to perform.
#' @param verbose Logical, if TRUE the transformation parameters are printed.
#' @param cofactor cofactor used for 'arcsinh' transformation, default 5, can be set to 150 for HDFC data.
#' @return transformed flow cytometry dataset
#'
#' @export
transform_data <- function(keep,
                           transformation,
                           original_data,
                           verbose,
                           cofactor = 5){

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
      temp <- original_data[,dataNum,drop=F] / cofactor
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
                           paramName, "; using default fluor logicle transformation, be carefull with this parameter!"))
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
      if (verbose == TRUE) {

        print(paste0(paramName, " w= ",w," t= ",t))

      }

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
#' @param transformation Transformation to perform. Select one of the following: \code{"auto_logi"} (autologicle, recommended for flow cytometry data), \code{"arcsinh"} (arcsinh transformation with cofactor 5), \code{"clr"} (centered-log-ratio) or \code{"none"} (no transformation).
#' @param remove_param Parameters to be removed from the \code{fcd}, "inTime" should be kept.
#' @param anno_table Path to the annotation table text file. The annotation table should contain one column with the file names of all .fcs or .csv files to read in and optionally additional columns with further sample information (e.g. "sample_id", "condition").
#' @param filename_col Name of the column of the \code{anno_table} containing the file name matching with the .fcs/.csv files.
#' @param seed A seed is set for reproducibility.
#' @param separator_anno Separator used in the annotation file, by default \code{separator_anno = ","}.
#' @param separator_fc_csv Separator used in the cytometry data .csv files, by default \code{separator_anno = ","}.
#' @param simple_names If TRUE only the channel description is used to name the column, if FALSE both channel name and description are pasted together.
#' @param truncate_max_range From FlowCore: logical type. Default is FALSE. can be optionally turned off to avoid truncating the extreme positive value to the instrument measurement range .i.e.'$PnR'.
#' @param emptyValue From FlowCore: boolean indicating whether or not we allow empty value for keyword values in TEXT segment. It affects how the double delimiters are treated. IF TRUE, The double delimiters are parsed as a pair of start and end single delimiter for an empty value. Otherwise, double delimiters are parsed one part of string as the keyword value. default is TRUE.
#' @param ignore.text.offset From FlowCore: whether to ignore the keyword values in TEXT segment when they don't agree with the HEADER. Default is FALSE, which throws the error when such discrepancy is found. User can turn it on to ignore TEXT segment when he is sure of the accuracy of HEADER so that the file still can be read.
#' @param verbose Default FALSE, if TRUE the at each file loaded something is printed in the screen.
#' @param cross_path_with_anno Defautl FALSE. If TRUE is the 'data_path' contains more files then the annotation table only the overlap will be loaded.
#' @param cofactor cofactor used for 'arcsinh' transformation, default 5, can be set to 150 for HDFC data.
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
                     transformation = NULL,
                     remove_param = NULL,
                     anno_table,
                     filename_col,
                     seed = 91,
                     separator_anno = ",",
                     separator_fc_csv = ",",
                     simple_names = TRUE,
                     truncate_max_range = FALSE,
                     emptyValue = TRUE,
                     ignore.text.offset = FALSE,
                     verbose = FALSE,
                     cross_path_with_anno = FALSE,
                     cofactor = 5) {

  # Set seed for reproducibility
  set.seed(seed)

  if (verbose) {

    print("Start reading the data")

  }

  ## Load the data
  data <- read_data(data_path = data_path,
                    max_cells = max_cell,
                    useCSV = useCSV,
                    separator = separator_fc_csv,
                    simple_names = simple_names,
                    truncate_max_range = truncate_max_range,
                    emptyValue = emptyValue,
                    ignore.text.offset = ignore.text.offset,
                    verbose = verbose,
                    cross_path_with_anno = cross_path_with_anno,
                    anno_table = anno_table,
                    separator_anno = separator_anno,
                    filename_col = filename_col)

  raw_data <- as.matrix(data$merged_df) # Take the dataframe with the intensity values

  ## Add InFile to remove_param

  remove_param <- c(remove_param, "InFile")

  ## Data Transformation
  keep <- colnames(raw_data)[!colnames(raw_data) %in% remove_param]

  raw_data <- raw_data[,which(colnames(raw_data) %in% keep)]

  ## Check if transformation parameter is provided
  if(!is.null(transformation)){
    ## Check if transformation is a valid value
    if (!transformation %in% c("clr", "arcsinh", "auto_logi", "none")) {
      stop(paste0(transformation, " is not a valid transformation method"))
    }
  }else{stop("transformation parameter needs to be specified to run this function")}

  if (verbose) {

    print("Start transforming the data")

  }



  trans_data <- transform_data(keep = keep,
                               transformation = transformation,
                               original_data = raw_data,
                               verbose = verbose,
                               cofactor = cofactor)

  ## Clean the dataframe
  df <- cbind(trans_data, expfcs_filename=data$merged_df[,"InFile"])
  df <- as.data.frame(df)
  df$expfcs_filename <- as.factor(df$expfcs_filename)
  df$expfcs_filename <- factor(df$expfcs_filename, labels = data$data_files)

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
  fcd[["extras"]][["prep_param"]] <- list(data_path = data_path,
                                          max_cell = max_cell,
                                          transformation = transformation,
                                          remove_param = remove_param,
                                          anno_table = anno_table,
                                          filename_col = filename_col,
                                          seed = seed,
                                          separator_anno = separator_anno,
                                          separator_fc_csv = separator_fc_csv,
                                          prep_function = "prep_fcd",
                                          version = packageDescription("cyCONDOR")$Version)

  class(fcd) <- "flow_cytometry_dataframe"

  return(fcd)

}


#' Read FlowJo workspace
#'
#' @title Read FlowJo Workspace
#' @description Loading and transforming the data to create a flow cytometry dataset from a Gate Set object for the analysis with the cyCONDOR workflow.
#' @param data_gs Gate Set object, e.g. created by using \code{\link[CytoML]{open_flowjo_xml}} and \code{\link[CytoML]{flowjo_to_gatingset}} from the CytoML package.
#' @param inverse.transform Logical: if the data should be reverse transformed or kept with FlowJo transformation (default = FALSE).
#' @param transformation If \code{inverse.transform = TRUE}, type of new transformation to perform. Select one of the following: \code{"auto_log"} (autologicle, recommended for flow cytometry data), \code{"arcsinh"} (arcsinh transformation), \code{"clr"} (centered-log-ratio) or \code{"none"} (no transformation).
#' @param remove_param Parameters to be removed from the \code{fcd}.
#' @param merge_anno Logical: If sample anno should be merged to the \code{fcd}.
#' @param anno_table If \code{merge_anno = TRUE}, path to the annotation table text file. The annotation table should contain one column with the file names of all .fcs or .csv files to read in and optionally additional columns with further sample information (e.g. "sample_id", "condition").
#' @param separator_anno Separator used in the annotation file, by default \code{separator_anno = ","}.
#' @param filename_col Name of the column of the \code{anno_table} containing the file name matching with the .fcs files.
#' @description \code{prep_fjw} takes a gate set object as input and returns a \code{fcd} with the FlowJo gating information saved in \code{fcd$anno$cell_anno}.
#' @return An object of class 'flow cytometry dataframe' (\code{fcd}) is returned.
#' @import flowWorkspace
#' @import Biobase
#' @import CytoML
#' @import stringr
#' @return read_flowjo_workspace
#'
#' @export
prep_fjw <- function(data_gs,
                     inverse.transform = FALSE,
                     transformation = NULL,
                     remove_param = NULL,
                     merge_anno = FALSE,
                     anno_table = NULL,
                     separator_anno = ",",
                     filename_col = NULL) {

  #get node list
  gate_list <- flowWorkspace::gs_get_pop_paths(data_gs, path = "auto")

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

    trans_data <- transform_data(keep = keep, transformation = transformation, original_data = raw_data, verbose = TRUE)

    ## Clean the dataframe
    df <- cbind(trans_data, data[, !colnames(data) %in% fs[[1]]@parameters$desc])
    df <- as.data.frame(df)
    colnames(df)[colnames(df) == "InFile"] <- "expfcs_filename"
    ## remove DIVA ID from file name
    filenames_clean <- str_split_fixed(filenames, pattern = ".fcs", n = 2)[, 1] %>%
      paste0(".fcs")
    df$expfcs_filename <- factor(df$expfcs_filename, labels = filenames_clean)

  } else {
    raw_data <- as.data.frame(data)
    ## filter out remove_param
    keep <- colnames(raw_data)[!colnames(raw_data) %in% remove_param]
    df <- raw_data[,which(colnames(raw_data) %in% keep)]
    colnames(df)[colnames(df) == "InFile"] <- "expfcs_filename"
    ## remove DIVA ID from file name
    filenames_clean <- str_split_fixed(filenames, pattern = ".fcs", n = 2)[, 1] %>%
      paste0(".fcs")
    df$expfcs_filename <- factor(df$expfcs_filename, labels = filenames_clean)


  }

  if (merge_anno == TRUE) {

    ## Now add the annotation (as csv file)
    anno <- read.delim(anno_table, sep = separator_anno)

    ## Check if file names in anno_table match those in df
    if(sum(anno[,filename_col] %in% df$expfcs_filename) < length(anno[,filename_col])){
      stop(paste0("The filenames in the annotation table ", anno_table, " do not match the file names present in your Gate Set object."))
    }

    df <- merge(df, anno, by.x = "expfcs_filename", by.y = filename_col)

  }

  ## Give the unique rownames
  rownames(df) <- paste(df$expfcs_filename, rownames(df), sep = "_")

  ## Prepare the final object
  fcd <- list()

  fcd[["expr"]][["orig"]] <- df[ ,colnames(df) %in% fs[[1]]@parameters$desc]
  fcd[["anno"]][["cell_anno"]] <- df[ ,!colnames(df) %in% fs[[1]]@parameters$desc]

  #save import parameters
  fcd[["extras"]][["prep_param"]] <- list(max_cell = 1000000000, #set to a high number to read in all events
                                          inverse.transform = inverse.transform,
                                          transformation = transformation,
                                          remove_param = remove_param,
                                          merge_anno = merge_anno,
                                          anno_table = anno_table,
                                          filename_col= filename_col,
                                          separator_anno = separator_anno,
                                          prep_function = "prep_fjw",
                                          version = packageDescription("cyCONDOR")$Version)


  class(fcd) <- "flow_cytometry_dataframe"

  return(fcd)

}


#' Save the results from cyCONDOR as FCS file(s)
#'
#' @title Save the results from cyCONDOR as FCS file(s)
#' @description Saves the expression data, annotation data and (optionally) results from the cyCONDOR analysis (dimensionality reduction, clustering) as one or more FCS file(s).
#' @param fcd flow cytometry data set
#' @param expr_slot expr_slot from which to take marker expression values, default is "orig".
#' @param reduction_method string specifying which dimensionality reduction method to use.
#' @param reduction_slot string specifying reduction name in reduction_method to use for visualization, e.g. "pca_orig".
#' @param cluster_slot string specifying which clustering slot to use to find variable specified in cluster_var.
#' @param cluster_var string specifying variable name in cluster_slot that identifies cell population labels to be used (e.g. clusters, metaclusters or predicted labels). Must be numeric in this function - factors are converted automatically.
#' @param split_by NULL or character string specifying a metadata variable (e.g. sample_ID) to split the data. If NULL, only one FCS file with all data is generated.
#' @param dir string specifying the directory where the FCS file(s) is/are saved. Current working directory by default.
#' @param filename string specifying the filename for the FCS file. If split_by is defined, an automatic suffix according to the defined variable is assigned to the FCS files.
#' @importFrom Biobase AnnotatedDataFrame
#' @import flowCore
#' @return Message that FCS file(s) are saved in the defined directory.
#'
#' @export
write_fcs <- function(fcd = condor,
                      expr_slot = "orig",
                      reduction_method = NULL,
                      reduction_slot = NULL,
                      cluster_slot = NULL,
                      cluster_var = NULL,
                      split_by = NULL,
                      dir = paste0(getwd(), "/"),
                      filename){

  # Print a working to let people know this is an experimental functio
  warning("This function is still experimental. Feel free to test if and report your experience, this function
  is currently writing cyCONDOR transformed value, future revision will export original values.")

  # Get the expression data
  xdata <- fcd[["expr"]][[expr_slot]]
  exp <- xdata

  # Add dimensionl reduction information
  if(!is.null(reduction_method)){
    rdata <- fcd[[reduction_method]][[reduction_slot]]
    exp <- cbind(exp, rdata)
  }

  # Add cluster labels (must be numeric!)
  if(!is.null(cluster_slot)){
    cluster <- fcd[["clustering"]][[cluster_slot]][[cluster_var]]
    exp <- cbind(exp, cluster)
  }

  # Ensure all columns in exp are numeric
  exp[] <- lapply(exp, function(x){
    if (is.numeric(x)) return(x)
    xnum <- as.numeric(as.character(x))
    if (any(is.na(xnum) & !is.na(x)))
      warning("Some values in column could not be converted to numeric.")
    return(xnum)
  })

  # Hanlde splitting
  if(!is.null(split_by)){
    anno <- fcd[["anno"]][["cell_anno"]]

    if(!split_by %in% colnames(anno)){
      stop("Column '", split_by, "' not found in the annotation data.")
    }

    split_vec <- anno[[split_by]]
    exp_mat <- as.matrix(exp)
    split_list <- split(seq_len(nrow(exp)), split_vec)

    # create a flowSet with flowFrames for each subset
    ffs <- lapply(names(split_list), function(group){
      idx <- split_list[[group]]
      exp_sub <- exp_mat[idx, , drop = FALSE]

      # Create AnnotatedDataFrame for parameters of each subset
      param_sub <- Biobase::AnnotatedDataFrame(data.frame(
        name = colnames(exp_sub),
        desc = paste("Channel", colnames(exp_sub), sep = "_"),
        range = apply(exp_sub, 2, function(x) ceiling(max(x))),
        minRange = apply(exp_sub, 2, min),
        maxRange = apply(exp_sub, 2, max)
      ))

      # flowFrames
      flowCore::flowFrame(exprs = exp_sub, parameters = param_sub)
    })

    names(ffs) <- names(split_list)

    # flowSet
    fs <- flowCore::flowSet(ffs)

    # Write each frame into its own fcs file
    for(sample_name in flowCore::sampleNames(fs)){
      out_file <- paste0(dir, filename, "_", make.names(sample_name), ".fcs")
      flowCore::write.FCS(fs[[sample_name]], file = out_file)
      message("File was saved to ", out_file)
    }
  }

  # No splitting, yields a single fcs file
  else{

    # Create AnnotatedDataFrame for parameters
    param <- Biobase::AnnotatedDataFrame(data.frame(
      name = colnames(exp),
      desc = paste("Channel", colnames(exp), sep = "_"),
      range = apply(exp, 2, function(x) ceiling(max(x))),
      minRange = apply(exp, 2, min),
      maxRange = apply(exp, 2, max)
    ))

    # Create flowFrame file
    ffs <- flowCore::flowFrame(exprs = as.matrix(exp),
                               parameters = param)

    # Write the fcs
    flowCore::write.FCS(ffs, file = paste0(dir, filename, ".fcs"))
    message("File was saved to ", dir, filename, ".fcs")
  }
}

#' Save the results from cyCONDOR as a Single Cell Experiment Class as RDS file
#'
#' @title Export Flow Cytometry Dataset as SCE.
#' @description This function creates an object of the class Single Cell Experiment (SCE) based on the data of the flow cytometry dataset (fcd).
#' @param fcd flow cytometry data set
#' @param expr_slot expr_slot from which to take marker expression values, default is "orig".
#' @param reduction_method string specifying which dimensionality reduction method to use.
#' @param reduction_slot string specifying reduction name in reduction_method to use for visualization, e.g. "pca_orig".
#' @param cluster_slot string specifying which clustering slot to use to find variable specified in cluster_var.
#' @param cluster_var string specifying variable name in cluster_slot that identifies cell population labels to be used (e.g. clusters, metaclusters or predicted labels). Must be numeric in this function - factors are converted automatically.
#' @param save boolean indicating if the sce should be saved as RDS file.
#' @param dir string specifying the directory where the RDS file is saved. Current working directory by default.
#' @param filename string specifying the filename for the RDS file.
#' @importFrom SingleCellExperiment SingleCellExperiment reducedDim
#' @importFrom SummarizedExperiment colData
#' @importFrom S4Vectors metadata
#' @returns The SCE  plus a short message if the RDS file was saved successfully.
#'
#' @export
export_sce <- function(fcd = condor,
                       expr_slot = "orig",
                       reduction_method = NULL,
                       reduction_slot = NULL,
                       cluster_slot = NULL,
                       cluster_var = NULL,
                       save = TRUE,
                       dir = paste0(getwd(), "/"),
                       filename = "") {

  message("This function is still experimental, if you run into any issue please let us know!")

  # Extracting the expression data
  exp <- as.matrix(t(fcd[["expr"]][[expr_slot]]))

  # Extracting the annotation data
  anno <- fcd[["anno"]][["cell_anno"]]
  names(anno)[names(anno) == "sample_ID"] <- "sample_id" # required for CATALYST - .check_sce() checks if 'sample_id' (case sensitive!) is present in the SCE.
  names(anno)[names(anno) == "group"] <- "condition" # also required naming for CATALYST

  # Check if the colnames of exp and the rownames of anno match
  if(!identical(colnames(exp), rownames(anno))){
    stop("Column names of the expression matrix and rownames of the annotation matrix do not match.")
  }

  # Create the SCE
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(assays = exp),
                                                    colData = anno)

  # Extract and add dimensional reduction data
  if(!is.null(reduction_method)){
    red <- fcd[[reduction_method]][[reduction_slot]]
    SingleCellExperiment::reducedDim(sce, reduction_method) <- red
  }

  # Extract and add clustering information
  if(!is.null(cluster_slot)){
    clust <- fcd[["clustering"]][[cluster_slot]][[cluster_var]]
    SummarizedExperiment::colData(sce)[[cluster_slot]] <- clust
    SummarizedExperiment::colData(sce)[["cluster_id"]] <- clust # also required for downstream application in CATALYST
  }

  # Extract and save experiment_info as metadata
  ei <- unique(as.data.frame(colData(sce)[, c("sample_id", "condition")]))
  rownames(ei) <- ei$sample_id  # Set rownames = sample IDs
  S4Vectors::metadata(sce)$experiment_info <- ei

  # Save the SCE
  ## as .rds for standard usage in R
  if(save == TRUE) {
    saveRDS(sce,
            file = paste0(dir, filename, ".rds"))
    message("SCE has been saved to ", dir, filename, ".rds")
  }

  return(sce)
}
