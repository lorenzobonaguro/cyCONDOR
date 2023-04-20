#' nfTransform
#'
#' @title nfTransform
#' @description Data transformation.
#' @param transTypeTable Table with the transformation parameters
#' @param dataA dataA
#' @param dataB dataB
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
#' @description Load the .fcs files into a dataframe
#' @param LoaderPATH Path to the .fcs files
#' @param ceil number of cells to subset
#' @param useCSV Logical, if input is .csv and not .fcs
#' @param separator Separato used the flow csv files (if loading from csv)
#' @import flowCore
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

#' plot_marker
#'
#' @title plot_marker
#' @description Function to visualize a selected parameter.
#' @param data This is the input dataframe for the visualization, this part of the code can still be improved but for the moment you have to cbind the dataframe with the informations you want to plot.
#' @param param Parameter to visualize in the plot, this can be either a continous variable or a categorical one, the function will react differently according.
#' @param order Logical if you want to order the dots in the plot, by expression for example. This can help to find small populations of positive cells.
#' @param title Title of the plot.
#' @param limX Limits of the x axes.
#' @param limY Limits of the y axes.
#' @param dim_red Dimensionality reduction to use for the visualization.
#' @param dot_size Size of the dots.
#' @param apha Alpha of the dots.
#' @param color_discrete Colors for discrete parameters.
#' @param color_gradient Colors for continous parameters.
#' @param remove_guide Logical, if you want to remove the guide.
#' @param facet_by_variable Logical if the plot should be splitted by the categorical variable used.
#' @import ggplot2
#' @import RColorBrewer
#' @import devtools
#' @import ggpubr
#' @import ggsci
#' @import ggrastr
#' @return plot marker or list of markers
#'
#' @export
plot_marker <- function(data,
                        param,
                        order = FALSE,
                        title = "adjust the title",
                        limX = NULL,
                        limY = NULL,
                        dim_red,
                        dot_size = 0.1,
                        apha = 0.2,
                        color_discrete = cluster_palette,
                        color_gradient = colors,
                        remove_guide = FALSE,
                        facet_by_variable = FALSE) {

  if(is.numeric(data[[param]]) == TRUE){

    if(order == TRUE){
      data <- data[order(data[[param]], decreasing = F), ]
    }

    colnames(data)[colnames(data) == param] <- "poi"

    if (dim_red == "UMAP") {

      p1 <- ggplot(data, aes(x = UMAP1, y = UMAP2, color = poi)) +
        geom_point(alpha = apha, size = dot_size) +
        theme_bw() +
        theme(aspect.ratio = 1, panel.grid = element_blank()) +
        ggtitle(title) + scale_color_gradientn(colours = color_gradient) + labs(color = param)

    }

    if (dim_red == "DM") {

      p1 <- ggplot(data, aes(x = DC_1, y = DC_2, color = poi)) +
        geom_point(alpha = apha, size = dot_size) +
        theme_bw() +
        theme(aspect.ratio = 1, panel.grid = element_blank()) +
        ggtitle(title) + scale_color_gradientn(colours = color_gradient) + labs(color = param)

    }

  }else {
    colnames(data)[colnames(data) == param] <- "poi"

    if (dim_red == "UMAP") {

      p1 <- ggplot(data, aes(x = UMAP1, y = UMAP2, color = poi)) +
        geom_point(alpha = apha, size = dot_size) +
        theme_bw() +
        theme(aspect.ratio = 1, panel.grid = element_blank()) +
        ggtitle(title) + scale_colour_manual(values = color_discrete) +
        guides(color = guide_legend(override.aes = list(size=5, alpha = 1))) + labs(color = param)

    }

    if (dim_red == "DM") {

      p1 <- ggplot(data, aes(x = DC_1, y = DC_2, color = poi)) +
        geom_point(alpha = apha, size = dot_size) +
        theme_bw() +
        theme(aspect.ratio = 1, panel.grid = element_blank()) +
        ggtitle(title) + scale_colour_manual(values = color_discrete) +
        guides(color = guide_legend(override.aes = list(size=5, alpha = 1))) + labs(color = param)

    }

    if (facet_by_variable == TRUE) {

      p1 <- p1 + facet_wrap(~poi)

    }


  }

  # Add some extra features on the plot
  if (!is.null(limX)) {

    p1 <- p1 + scale_x_continuous(limits = limX)

  }

  if (!is.null(limY)) {

    p1 <- p1 + scale_y_continuous(limits = limY)

  }

  if (remove_guide == TRUE) {

    p1 <- p1 + theme(legend.position = "none")

  }


  return(p1)

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

  for (level_1 in names(fcdataset)) {

    int_cont <- fcdataset[[level_1]]

    int_collector <- list()

    for (level_2 in names(int_cont)) {

      int_collector[[level_2]] <- int_cont[[level_2]][cell_ids, ]
    }

    new_fcd[[level_1]] <- int_collector

  }

  return(new_fcd)

}

#' check_IDs
#'
#' @title check_IDs
#' @description CHeck the integrity and the correctness of a flow cytometry dataset
#' @param fcdataset Flow cytometry dataset to be checked.
#' @return check the integrity and the correctnes of a flow cytometry dataset
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

#' scaleColors
#'
#' @title scaleColors
#' @description Defines the color coding for an heatmap.
#' @param data data to use
#' @param maxvalue value at which the color is fully red / blue
#' @return Defines the color coding for an heatmap.
#'
#' @export
scaleColors <- function(data = input_scale,
                        maxvalue = NULL
){
  if(is.null(maxvalue)){
    maxvalue <- floor(min(abs(min(data)), max(data)))
  }
  if(max(data) > abs(min(data))){
    if(ceiling(max(data)) == maxvalue){
      myBreaks <- c(floor(-max(data)), seq(-maxvalue+0.2, maxvalue-0.2, 0.2),  ceiling(max(data)))
    } else{
      myBreaks <- c(floor(-max(data)), seq(-maxvalue, maxvalue, 0.2),  ceiling(max(data)))
    }
    paletteLength <- length(myBreaks)
    myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
  } else {
    if(-floor(min(data)) == maxvalue){
      myBreaks <- c(floor(min(data)), seq(-maxvalue+0.2, maxvalue-0.2, 0.2),  ceiling(min(data)))
    } else{
      myBreaks <- c(floor(min(data)), seq(-maxvalue, maxvalue, 0.2),  ceiling(abs(min(data))))
    }
    paletteLength <- length(myBreaks)
    myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
  }
  return(list(breaks = myBreaks, color = myColor))
}

#' confusionMatrix
#'
#' @title confusionMatrix
#' @description Calculate a confusion matrix
#' @param i i
#' @param j j
#' @importFrom Matrix sparseMatrix
#' @import pheatmap
#' @return Create a confusion matrix
#'
#' @export
confusionMatrix <- function (i = NULL, j = NULL)
{
  ui <- unique(i)
  uj <- unique(j)
  m <- Matrix::sparseMatrix(i = match(i, ui), j = match(j,
                                                        uj), x = rep(1, length(i)), dims = c(length(ui), length(uj)))
  rownames(m) <- ui
  colnames(m) <- uj
  m
}

#' HM_markers
#'
#' @title HM_markers
#' @description Visualize the heatmap of marker expression groupped by a variable of interest.
#' @param input cbind of the expression table to be used the the groupping (eg. clustering).
#' @param group column name to be used for the groupping (eg. "Phenograph" or "group")
#' @param maxvalue -Max scaled expression to be used for the color coding.
#' @param size Size of the cells of the heatmap.
#' @param title Title for the plot.
#' @param exclusion Marker to exclude from the visualization.
#' @import reshape2
#' @return HM_markers
#'
#' @export
HM_markers <- function(input,
                       group,
                       maxvalue = NULL,
                       size = 10,
                       title = "Choose a good title for the plot",
                       exclusion = c("SSC-A", "FSC-A"),
                       cluster_rows = FALSE,
                       cluster_cols = FALSE) {

  tmp <- suppressMessages(melt(input))

  tmp <- summarySE(data = tmp, measurevar = "value", groupvars = c(group, "variable"))

  colnames(tmp)[colnames(tmp) == group] <- "group"

  tmp <- dcast(tmp[, c(1,2,4)], variable ~ group)

  rownames(tmp) <- tmp$variable
  tmp$variable <- NULL
  tmp <- tmp[!rownames(tmp) %in% exclusion,]

  tmp <- t(scale(t(tmp)))

  pheatmap::pheatmap(tmp, scale = "none",
                     breaks = scaleColors(data = tmp, maxvalue = maxvalue)[["breaks"]],
                     color = scaleColors(data = tmp, maxvalue = maxvalue)[["color"]],
                     main = title, cellwidth = size, cellheight = size,
                     cluster_rows = cluster_row, cluster_cols = cluster_cols)

}

#' confusion_HM
#'
#' @title confusion_HM
#' @description Visualize the confusion heatmap.
#' @param variables Variable to by used to calculate the confusion.
#' @param group Groupping to calculate the relative contribution to the variable.
#' @param size Size of the heatmap blocks.
#' @param title Title of the plot.
#' @return confusion_HM
#'
#' @export
confusion_HM <- function(variables, group, size = 15, title = "Choose a good title for the plot") {

  # quantify cells of each sample per cluster
  cells_cluster <- confusionMatrix(paste0(variables),
                                   paste0(group))

  cells_cluster <- cells_cluster[order(factor(rownames(cells_cluster),levels=c(0:nrow(cells_cluster)))),]

  cells_cluster <- cells_cluster[, order(colnames(cells_cluster))]

  cells_cluster <- as.matrix(cells_cluster)

  # normalize to 1000 cells per sample
  tmp <- round(t(t(cells_cluster)/colSums(cells_cluster))*1000,3)
  # calculate percentage of cells from sample per cluster
  scaled_cM <- round((tmp / Matrix::rowSums(tmp))*100,2)

  pheatmap::pheatmap(
    mat = t(scaled_cM),
    border_color = "black",display_numbers = TRUE,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    cellwidth = size,
    cellheight = size,
    # breaks = c(seq(0, maxvalue, length.out = 100), 100),
    #color = c(colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100), "#4575B4"),
    main = title)

}

#' boxplot_and_stats
#'
#' @title boxplot_and_stats
#' @description Calculate boxplot and summary statistics for a wanted comparison.
#' @param annotation Sample annotation to be used for the plot.
#' @param sample_var Column name containing the sample IDs.
#' @param group_var Column name defining the groupping for plotting.
#' @param variable Variable used to stratify the plotting.
#' @param numeric Logical if the groupping is numeric.
#' @param test.type Test to be performed. (see need some development here)
#' @return boxplot_and_stats
#'
#' @export
boxplot_and_stats <- function(annotation,
                              sample_var,
                              group_var,
                              variable,
                              numeric = TRUE,
                              test.type = "t.test") {

  # Container
  container <- list()

  # quantify cells of each sample per cluster
  tmp <- confusionMatrix(paste0(annotation[[sample_var]]),
                         paste0(variable))

  tmp <- as.matrix(tmp)

  if (numeric == TRUE) {

    tmp <- tmp[order(rownames(tmp)),order(as.numeric(colnames(tmp)))]

  }

  if (numeric == FALSE) {

    tmp <- tmp[order(rownames(tmp)),order(colnames(tmp))]

  }

  # Give a more complete naming to the columns
  colnames(tmp) <- paste("Cluster_", colnames(tmp), sep = "")

  # calculate the percentage
  tmp <- tmp/rowSums(tmp)*100

  tmp <- as.data.frame(tmp)

  # Merging with the annotation
  tmp[[sample_var]] <- rownames(tmp)

  tmp <- merge(tmp, unique(annotation), by = sample_var)

  container[["per_data"]] <- tmp

  # Plotting
  colnames(tmp)[colnames(tmp) == group_var] <- "groups"

  for (variable in colnames(tmp)[grep("Cluster_", colnames(tmp))]) {

    data <- tmp

    colnames(data)[colnames(data) == variable] <- "poi"

    p1 <- ggplot(data, aes(x = groups, y = poi, fill = groups)) +
      geom_boxplot(outlier.colour = NA) +
      geom_point(shape = 21, fill = "white", color = "black", size = 4, position = position_dodge(0.4)) +
      theme_bw() + theme(aspect.ratio = 2, panel.grid = element_blank()) +
      scale_fill_manual(values = cluster_palette) +
      ggtitle(variable) +
      expand_limits(y = 0) +
      xlab("") +
      stat_compare_means(method = test.type, paired = FALSE, label = "p.signif", label.x.npc = 0.5)

    container[["plot"]][[variable]] <- p1

    rm(p1, data)

  }

  tmp2 <- suppressMessages(melt(tmp[, c(colnames(tmp)[grep("Cluster_", colnames(tmp))], "groups")]))

  stat <- compare_means(formula = value~groups,
                        data = tmp2,
                        paired = TRUE,
                        method = test.type,
                        group.by = "variable",
                        p.adjust.method = "none")

  container[["stats"]] <- stat



  return(container)

}

#' barplot_frequency
#'
#' @title barplot_frequency
#' @description The function get_expr_huva provides the expression table resulting from reference huva experiment data.frame.
#' @param x_axes Groupping of the x axes.
#' @param colour Stratification to use on the stacked barplot.
#' @param color_palette color_palette
#' @param title Title for the plot.
#' @param legend_title Title for the legend.
#' @return barplot_frequency
#'
#' @export
barplot_frequency <- function(x_axes,
                              colour,
                              color_palette = cluster_palette,
                              title = "Choose an appropriate title",
                              legend_title) {

  tmp <- table(colour,
               x_axes)

  tmp <- t(t(tmp)/colSums(tmp))*100

  tmp <- as.data.frame(tmp)

  p1 <- ggplot(tmp, aes(x = x_axes, y = Freq, fill = colour)) +
    geom_bar(stat="identity", color = "black", size = 0.3) +
    scale_fill_manual(values = color_palette) +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    labs(y = "Frequency (%)", x = "", fill = legend_title) +
    ggtitle(title)

  return(p1)

}

#' prep_fcd
#'
#' @title prep_fcd
#' @description Wrapping function to prepare a flow cytometry dataset
#' @param FCSpath Folder where the .fcs files are stored.
#' @param ceil Number of cells to use for each file (set to a high number if you want to use all available events)
#' @param useCSV Flag if the input are .csv files and not .fcs (experimental).
#' @param transformation Transformation to perform.
#' @param remove_param Parameters to remove from the trasfomration, "inTime" should be kept.
#' @param anno_table path to the annotation table file.
#' @param filename_col Name of the column containing the filename matching with the .fcs files.
#' @param seed seed to be used for the randomization of the events.
#' @param separator_anno separator used in the annotation file
#' @param separator_fc_csv separator used in the fc csv files
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

#' harmonize_intensities
#'
#' @title harmonize_intensities
#' @description Harmonize the expression values.
#' @param fcd flow cytometry dataset
#' @param batch vector of column names to use for correcting the data
#' @param seed Seed used for the randomization steps
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

#' runPCA
#'
#' @title runPCA
#' @description Run a Principal Component Analysis.
#' @param fcd flow cytometry dataset
#' @param data_slot name of the data slot to use to calculate the PCA, original data (orig) or harmonized data (norm)
#' @param seed Seed used for the randomization steps.
#' @param prefix Prefix for the output.
#' @return runPCA
#'
#' @export
runPCA <- function(fcd, data_slot = "orig", seed, prefix = NULL) {

  set.seed(seed)

  if (is.null(prefix)) {

    fcd[["pca"]][[data_slot]] <- prcomp(fcd$expr[[data_slot]])$x

  } else {

    fcd[["pca"]][[paste(prefix, data_slot, sep = "_")]] <- prcomp(fcd$expr[[data_slot]])$x

  }

  return(fcd)

}

#' harmonize_PCA
#'
#' @title harmonize_PCA
#' @description Harmonize the Principal Component Analysis.
#' @param fcd flow cytometry dataset
#' @param data_slot name of the PCA data slot to use to harmonize. If no prefix was added the, *orig*.
#' @param batch vector of column names to use for correcting the data
#' @param seed Seed used for the randomization steps
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

#' runUMAP
#'
#' @title runUMAP
#' @description Run a UMAP dimensionality reduction.
#' @param fcd flow cytometry dataset
#' @param input_type data to use for the calculation of the UMAP, e.g. "expr" or "pca"
#' @param data_slot name of the PCA data slot to use to harmonize. If no prefix was added the, *orig*.
#' @param n_neighbors n_neighbors
#' @param n_components n_components
#' @param min_dist min_dist
#' @param metric metric
#' @param seed Seed used for the randomization steps
#' @param prefix Prefix for the output.
#' @param n_threads Number of threads to be used in the UMAP calculation.
#' @import umap
#' @import Rtsne
#' @return runUMAP
#'
#' @export
runUMAP <- function (fcd,
                     input_type,
                     data_slot,
                     n_neighbors = 15,
                     n_components = 2,
                     min_dist = 0.2,
                     metric = "euclidean",
                     seed,
                     prefix = NULL,
                     n_threads = 32) {

  set.seed(seed)

  umapMat <- uwot::umap(X = fcd[[input_type]][[data_slot]],
                        n_neighbors = n_neighbors,
                        n_components = n_components,
                        min_dist = min_dist,
                        metric = metric,
                        n_threads = n_threads)

  colnames(umapMat) <- c("UMAP1", "UMAP2")

  if (is.null(prefix)) {

    fcd[["umap"]][[paste(input_type, data_slot, sep = "_")]] <- umapMat

  }else {

    fcd[["umap"]][[paste(prefix, input_type, data_slot, sep = "_")]] <- umapMat
  }

  return(fcd)
}

#' runPhenograph
#'
#' @title runPhenograph
#' @description Run Phenograph clustering.
#' @param fcd flow cytometry dataset.
#' @param input_type data to use for the calculation of the UMAP, e.g. "pca" (suggested option).
#' @param data_slot name of the PCA data slot to use to harmonize. If no prefix was added the, *orig*.
#' @param k K value used for clustering.
#' @param seed Seed used for the randomization steps.
#' @param prefix Prefix for the output.
#' @import Rphenograph
#' @importFrom igraph membership
#' @return runPhenograph
#'
#' @export
runPhenograph <- function(fcd, input_type, data_slot, k, seed, prefix = NULL) {

  set.seed(seed)

  Rphenograph_out <- Rphenoannoy::Rphenoannoy(fcd[[input_type]][[data_slot]], k = k)
  Rphenograph_out <- as.matrix(membership(Rphenograph_out[[2]]))
  Rphenograph_out <- as.data.frame(matrix(ncol=1,data=Rphenograph_out,dimnames=list(rownames(fcd$expr$orig),"Phenograph")))

  Rphenograph_out$Phenograph <- as.factor(Rphenograph_out$Phenograph)
  Rphenograph_out$Description <- paste(input_type, "_",data_slot, "_k", k, sep = "")

  if (is.null(prefix)) {

    fcd[["clustering"]][[paste("Phenograph_", input_type, "_",data_slot, "_k", k, sep = "")]] <- Rphenograph_out

  } else {

    fcd[["clustering"]][[paste("Phenograph_", prefix, "_",input_type, "_",data_slot, "_k", k, sep = "")]] <- Rphenograph_out

  }

  return(fcd)

}

#' runDM
#'
#' @title runDM
#' @description Run Diffusion Map dimensionality reduction.
#' @param fcd flow cytometry dataset
#' @param input_type data to use for the calculation of the UMAP, e.g. "pca" (suggested option)
#' @param data_slot name of the PCA data slot to use to harmonize. If no prefix was added the, *orig*.
#' @param k K used for the analysis.
#' @param seed Seed used for the randomization steps
#' @param prefix Prefix of the output.
#' @import destiny
#' @import SingleCellExperiment
#' @import slingshot
#' @return runDM
#'
#' @export
runDM <- function(fcd, input_type, data_slot, k = 10,seed, prefix = NULL) {

  set.seed(91)

  dm <- DiffusionMap(fcd[[input_type]][[data_slot]],
                     vars = NULL,
                     k=k, suppress_dpt = TRUE, verbose=TRUE, n_pcs = NA)

  dm <- cbind(dm$DC1, dm$DC2, dm$DC3)
  colnames(dm) <- c("DC_1", "DC_2", "DC_3")

  if (is.null(prefix)) {

    fcd[["diffmap"]][[paste(input_type, data_slot, sep = "_")]] <- dm

  } else {

    fcd[["diffmap"]][[paste(prefix, input_type, data_slot, sep = "_")]] <- dm

  }

  return(fcd)

}

#' merge_condor
#'
#' @title merge_condor
#' @description Merges two condor objects.
#' @param data1 Dataset 1 to merge
#' @param data2 Dataset 2 to merge
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

#' metaclustering
#'
#' @title metaclustering
#' @description Assignment of a metacluster.
#' @param fcd flow cytometry dataset.
#' @param clustering Name of the clustering to match for the metaclustering.
#' @param name_col Column containing the original cluster
#' @param name_out Name of the output column
#' @param metaclusters Vector of the new clusters names, this should be of the same length of the levels of the original clustering.
#' @return metaclustering
#'
#' @export
metaclustering <- function(fcd,
                           clustering,
                           name_col,
                           name_out = "metacluster",
                           metaclusters) {

  ### Assign the metacluster
  #### To make it easier when the number of cluster is high
  tmp_metaclusters <- data.frame(cluster = levels(fcd$clustering[[clustering]][[name_col]]),
                                 metacluster = c(metaclusters))

  print(tmp_metaclusters)

  fcd$clustering[[clustering]][[name_out]] <- factor(fcd$clustering[[clustering]][[name_col]],
                                                     levels = levels(fcd$clustering[[clustering]][[name_col]]),
                                                     labels = tmp_metaclusters$metacluster)

  return(fcd)

}

#' runFlowSOM
#'
#' @title runFlowSOM
#' @description Run FlowSOM based clustering.
#' @param fcd flow cytometry dataset.
#' @param input_type data to use for the calculation of the UMAP, e.g. "expr" or "pca".
#' @param data_slot name of the PCA data slot to use to harmonize. If no prefix was added the, *orig*.
#' @param num_clusters number of final clusters.
#' @param seed Seed used for the randomization steps.
#' @param prefix Prefic for the output.
#' @return metaclustering
#'
#' @export
runFlowSOM <- function(fcd, input_type, data_slot, num_clusters, seed, prefix = NULL) {

  set.seed(seed)

  out <- FlowSOM::ReadInput(fcd[[input_type]][[data_slot]], transform = FALSE, scale = FALSE)

  out <- FlowSOM::BuildSOM(out, colsToUse = 1:(ncol(fcd[[input_type]][[data_slot]])))

  out <- FlowSOM::BuildMST(out)

  labels_pre <- out$map$mapping[, 1]

  out <- FlowSOM::metaClustering_consensus(out$map$codes, k = num_clusters, seed = seed)

  labels <- out[labels_pre]

  labels <- as.factor(labels)
  Description <- paste(input_type, "_",data_slot, "_k", num_clusters, sep = "")

  flSOM <- data.frame("FlowSOM" = labels, "Description" = Description)

  rownames(flSOM) <- rownames(fcd$expr$orig)

  if (is.null(prefix)) {

    fcd[["clustering"]][[paste("FlowSOM_", input_type, "_",data_slot, "_k", num_clusters, sep = "")]] <- flSOM

  } else {

    fcd[["clustering"]][[paste("FlowSOM_", prefix, "_",input_type, "_",data_slot, "_k", num_clusters, sep = "")]] <- flSOM

  }

  return(fcd)

}
