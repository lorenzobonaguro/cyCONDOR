#' nfTransform
#'
#' @title nfTransform
#' @description Data transformation, this function run within the prep_fcd wrapper.
#' @param transTypeTable Table with the transformation parameters.
#' @param dataA dataA.
#' @param dataB dataB, same as dataA.
#' @return transformed flow cytometry dataset
#'
#' @export
nfTransform <- function(transTypeTable, dataA, dataB){
  dataA1<-dataA
  dataB1<-dataB
  CyTOFlgcl <- logicleTransform(w=0.25, t=16409, m=4.5, a=0)
  ilgcl <- inverseLogicleTransform(trans = CyTOFlgcl)
  Fluorlgcl <- logicleTransform(w=0.1, t=500000, m=4.5, a=0)


  #perform transform
  for(paramName in as.character(transTypeTable[,1]))
  {
    ttParamNum <- which(as.character(transTypeTable[,1])==paramName)
    ttParamType <- as.character(transTypeTable[ttParamNum,2])
    if(ttParamType == "y" || ttParamType == "c"){
      dataNum <- which(colnames(dataA)==paramName)
      temp <- apply(dataA[,dataNum,drop=F],2, CyTOFlgcl)
      dataA1[,dataNum] <- temp
      dataB1[,dataNum] <- temp
      #print(paste(paramName," CyTOFlgcl"))
    }

    if(ttParamType == "f" ){
      dataNum <- which(colnames(dataA)==paramName)
      temp <- apply(dataA[,dataNum,drop=F],2, Fluorlgcl)
      dataA1[,dataNum] <- temp
      dataB1[,dataNum] <- temp
    }
    if(ttParamType == "l" ){
      dataNum <- which(colnames(dataA)==paramName)
      temp <- (dataA[,dataNum]/max(dataA[,dataNum]))*4.5
      dataA1[,dataNum] <- temp
      dataB1[,dataNum] <- temp
    }
    if(ttParamType == "n" ){
      dataNum <- which(colnames(dataA)==paramName)
      temp <- ((dataA[,dataNum]-min(dataA[,dataNum]))/(max(dataA[,dataNum])-min(dataA[,dataNum])))*4.5
      dataA1[,dataNum] <- temp
      dataB1[,dataNum] <- temp

    }
    if(ttParamType == "a"){
      q<-0.05
      m<-4.5
      d <- dataA[,paramName]
      w <- 0
      t <- max(d)
      nd <- d[d < 0]
      nThres <- quantile(nd, 0.25) - 1.5 * IQR(nd)
      nd <- nd[nd >= nThres]
      #transId <- paste(p, "autolgclTransform", sep = "_")
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
      dataNum <- which(colnames(dataA)==paramName)
      temp <- apply(dataA[,dataNum,drop=F],2, templgcl)
      dataA1[,dataNum] <- temp
      dataB1[,dataNum] <- temp
      print(paste0(paramName, " w= ",w," t= ",t))
      #hist(temp, main=paramName,breaks=100)
    }

  }
  return(list(dataA1=dataA1, dataB1=dataB1))
}

#' prepFcsFolderData
#'
#' @title prepFcsFolderData
#' @description Load .fcs or .csv files into a dataframe and prepare the condor object.
#' @param LoaderPATH Path to the .fcs files.
#' @param ceil number of cells to subset.
#' @param useCSV Logical, if input is .csv and not .fcs.
#' @param separator Separator used the flow csv files (if loading from csv).
#' @import flowCore
#' @import reshape2
#' @import dplyr
#' @import cowplot
#' @return load flow cytometry dataset
#'
#' @export
prepFcsFolderData <- function(LoaderPATH, ceil, useCSV, separator){
  if(!useCSV){



    FcsFileNames <- list.files(path = LoaderPATH, pattern = ".fcs")
    fs = list()
    for(FileNum in 1:length(FcsFileNames)){
      fs[[FileNum]] <- read.FCS(paste0(LoaderPATH,"/",FcsFileNames[FileNum]),transformation =FALSE,ignore.text.offset=T,truncate_max_range=FALSE,emptyValue = F)
    }


    #fs <- read.flowSet(path = LoaderPATH,transformation=FALSE,ignore.text.offset=T,truncate_max_range=FALSE,emptyValue = F)
    #FcsFileNames <- rownames(keyword(fs, "FILENAME"))
    filenames <- FcsFileNames
    NumBC <- length(fs)
    FFdata <- NULL
    OrigNames <-fs[[1]]@parameters$name
    for (FFs in 1:NumBC){
      FFt <- exprs(fs[[FFs]])
      ## Downsample
      if (nrow(FFt)<=ceil)
        FFa <- FFt
      else
        FFa <- FFt[sample(nrow(FFt),ceil,replace=F),]
      #Fixup column names
      colnames(FFa) <- fs[[FFs]]@parameters$desc
      empties <- which(is.na(colnames(FFa)) | colnames(FFa)== " ")
      colnames(FFa)[empties] <- fs[[FFs]]@parameters$name[empties]
      fs[[FFs]]@parameters$desc <- colnames(FFa)

      #fs[[FFs]]@parameters$name <- colnames(FFa)
      #Add file label
      FFa <- cbind(FFa,rep(FFs,dim(FFa)[1]))
      colnames(FFa)[dim(FFa)[2]] <- "InFile"
      #Concatenate
      FFdata <- rbind(FFdata,FFa)
    }} else {  #load csv files instead of fcs

      csvfilenames <- list.files(path = LoaderPATH, pattern="*.csv")
      FcsFileNames <- csvfilenames
      csvdata <- lapply(paste0(LoaderPATH,"//",csvfilenames),function(x) read.delim(x, check.names = F, sep = separator_fc_csv))
      NumBC <- length(csvdata)
      FFdata<-NULL
      for (FFs in 1:NumBC){
        FFt <- csvdata[[FFs]]
        ## Downsample
        if (nrow(FFt)<=ceil)
          FFa <- FFt
        else
          FFa <- FFt[sample(nrow(FFt),ceil,replace=F),]

        FFa <- cbind(FFa,rep(FFs,dim(FFa)[1]))
        colnames(FFa)[dim(FFa)[2]] <- "InFile"
        #Concatenate
        FFdata <- rbind(FFdata,FFa)
      }
    }
  if(!useCSV){

    return(list(FFdata=FFdata, OrigNames=OrigNames, forNewFF=fs[[1]], NumBC=NumBC, FcsFileNames = FcsFileNames))

  }else{

    return(list(FFdata=FFdata, NumBC=NumBC, FcsFileNames = FcsFileNames))
  }
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
                     transformation = 'a',
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
    keeptable <- data.frame(Param = fs[[1]]@parameters$desc)
    keeptable$Trans <- transformation
    keeptable <- keeptable[!keeptable$Param %in% remove_param, ]

    data1 <- FFdata[,which(colnames(FFdata) %in% keeptable[,1])]

    nfTransOut <- nfTransform(keeptable, data1, data1)

    data1 <- nfTransOut$dataA1

    ## Clean the dataframe
    df <- cbind(data1, FFdata[, !colnames(FFdata) %in% fs[[1]]@parameters$desc])
    df <- as.data.frame(df)
    colnames(df)[colnames(df) == "InFile"] <- "expfcs_filename"
    df$expfcs_filename <- as.factor(df$expfcs_filename)
    df$expfcs_filename <- factor(df$expfcs_filename, labels = filenames)

  } else {

    df <- FFdata
    colnames(df)[colnames(df) == "InFile"] <- "expfcs_filename"
    df$expfcs_filename <- as.factor(df$expfcs_filename)
    df$expfcs_filename <- factor(df$expfcs_filename, labels = filenames)

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

#' prep_fcd
#'
#' @title prep_fcd
#' @description Wrapping function to prepare a flow cytometry dataset
#' @param FCSpath Folder where the .fcs files are stored.
#' @param ceil Number of cells to use for each file (set to a high number if you want to use all available events).
#' @param useCSV Flag if the input are .csv files and not .fcs (experimental).
#' @param transformation Transformation to perform.
#' @param remove_param Parameters to remove from the transformation, "inTime" should be kept.
#' @param anno_table path to the annotation table file.
#' @param filename_col Name of the column containing the filename matching with the .fcs files.
#' @param seed seed to be used for the randomization of the events.
#' @param separator_anno separator used in the annotation file.
#' @param separator_fc_csv separator used in the fc csv files.
#' @import readr
#' @import readxl
#' @import stringr
#' @import Rmisc
#' @return prep_fcd
#'
#' @export
prep_fcd <- function(FCSpath,
                     ceil,
                     useCSV = FALSE,
                     transformation,
                     remove_param,
                     anno_table,
                     filename_col,
                     seed,
                     separator_anno = ",",
                     separator_fc_csv = ",") {

  # Set seed for reproducibility
  set.seed(seed)

  ## Load the data
  prepData <- prepFcsFolderData(LoaderPATH = FCSpath, ceil = ceil, useCSV = useCSV, separator = separator_fc_csv)

  FFdata<- as.matrix(prepData$FFdata) # Take the dataframe with the intensity values

  ## Data Transformation
  keeptable <- data.frame(Param = colnames(FFdata))
  keeptable$Trans <- transformation
  keeptable <- keeptable[!keeptable$Param %in% remove_param, ]

  data1 <- FFdata[,which(colnames(FFdata) %in% keeptable[,1])]

  nfTransOut <- nfTransform(keeptable, data1, data1)

  data1 <- nfTransOut$dataA1

  ## Clean the dataframe
  df <- cbind(data1, expfcs_filename=FFdata[,"InFile"])
  df <- as.data.frame(df)
  df$expfcs_filename <- as.factor(df$expfcs_filename)
  df$expfcs_filename <- factor(df$expfcs_filename, labels = prepData$FcsFileNames)

  ## Now add the annotation (as csv file)
  anno <- read.delim(anno_table, sep = separator_anno)

  df <- merge(df, anno, by.x = "expfcs_filename", by.y = filename_col)

  ## Give the unique rownames
  rownames(df) <- paste(df$expfcs_filename, rownames(df), sep = "_")

  ## Prepare the final object
  fcd <- list()

  fcd[["expr"]][["orig"]] <- df[ ,colnames(df) %in% keeptable$Param]
  fcd[["anno"]][["cell_anno"]] <- df[ ,!colnames(df) %in% keeptable$Param]

  class(fcd) <- "flow_cytometry_dataframe"

  return(fcd)

}
