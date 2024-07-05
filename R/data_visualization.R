#' getTable
#'
#' @title Get summary values
#' @description
#' \code{getTable()} can be used to generate frequently used parameters on cell populations defined by clustering or prediction, while considering a meta variable for grouping. It can produce cell numbers (counts), cell population frequencies as well as median or mean marker expression for each group to cell population combination.
#' @param fcd flow cytometry data set, that has been subjected to clustering or cell type label prediction with cyCONDOR
#' @param expr_slot expr_slot from which to take marker expression values, default is "orig".
#' Corrected input data should be handled cautiously.
#' @param output_type type of parameter that should be reported in table. One of the following option needs to be selected:
#' \itemize{
#'  \item{"counts"} : gives cell numbers per group_var and cell population
#'  \item{"frequency"} : returns proportion of each cell population for each level in group_var (default)
#'  \item{"median"}: calculates median expression for each group_var and cell population combination for each available feature in expression matrix
#'  \item{"mean"}: calculates mean expression for each group_var and cell population combination for each available feature in expression matrix
#' }
#' @param cluster_slot string specifying which clustering slot to use to find variable specified in cluster_var
#' @param cluster_var string specifying variable name in cluster_slot that identifies cell population labels to be used (e.g. clusters, metaclusters or predicted labels).
#' @param group_var string indicating variable name in cell_anno that should be used to group the output, e.g. group or sample ID.
#' @param numeric logical, if TRUE numeric levels in cluster_var are ordered in increasing order and "Cluster_" is pasted before number, if FALSE alphabetical ordering is applied.
#' @import dplyr
#' @returns
#' \code{getTable()} returns a data frame with parameters in columns and observations in rows. In case of output_type of "counts" or "frequency", counts and frequencies for each cell population (columns) are reported in one row for each level in group_var. Given an output_type of "mean" or "median", aggregated expression for each feature (columns) is reported for each group_var and cell population combination (cluster).
#'
#'@export
getTable <- function(fcd,
                     output_type = "frequency", #alternative "counts","mean","median"
                     expr_slot = "orig",
                     cluster_slot,
                     cluster_var,
                     group_var,
                     numeric = F) {


  if(!output_type %in% c("frequency","counts","mean","median")){
    stop('output_type must be set to "counts","frequency","mean" or "median".')

  }else if(output_type %in% c("counts","frequency")){

    #### check slots, cell IDs and variables
    checkInput(fcd = fcd,
               check_cluster_slot = T,
               check_cell_anno = T,
               cluster_slot = cluster_slot,
               cluster_var = cluster_var,
               group_var = group_var)


    #### prepare data
    data <- data.frame(cellID = rownames(fcd$clustering[[cluster_slot]]),
                       cluster = fcd$clustering[[cluster_slot]][[cluster_var]],
                       group_var = fcd$anno$cell_anno[[group_var]])


    #### prepare output
    counts <- cyCONDOR::confusionMatrix(paste0(data$group_var), paste0(data$cluster))
    counts <- as.matrix(counts)
    if (numeric == TRUE) {
      counts <- counts[order(rownames(counts)), order(as.numeric(colnames(counts)))]
      colnames(counts) <- paste0("Cluster_", colnames(counts))
    } else if (numeric == FALSE) {
      counts <- counts[order(rownames(counts)), order(colnames(counts))]
    } else {
      counts <- counts[order(rownames(counts)), order(colnames(counts))]
    }

    if(output_type == "counts"){
      counts <- as.data.frame(counts)
      counts$group_var <- rownames(counts)
      counts <- counts[,c("group_var", c(colnames(counts)[!colnames(counts) == "group_var"]))]
      return(counts)

    }else{
      frequency <- counts/rowSums(counts) * 100
      frequency <- as.data.frame(frequency)
      frequency$group_var <- rownames(frequency)
      frequency <- frequency[,c("group_var",colnames(frequency)[!colnames(frequency)=="group_var"])]
      return(frequency)

    }

  }else{
    checkInput(fcd = fcd,
               check_expr_slot = T,
               check_cluster_slot = T,
               check_cell_anno = T,
               expr_slot = expr_slot,
               cluster_slot = cluster_slot,
               cluster_var = cluster_var,
               group_var = group_var)

    if(is.null(expr_slot)){
      stop('To calculate mean or median expression, the expr_slot from which to take the expression values needs to be specified.')
    }

    #### prepare data
    data <- fcd$expr[[expr_slot]]
    data$cluster <- fcd$clustering[[cluster_slot]][[cluster_var]]
    data$group_var <- fcd$anno$cell_anno[[group_var]]


    ####calculate means per group x cluster combination
    if(output_type == "mean"){
      data_stats <- data %>% dplyr::group_by(group_var, cluster, .drop=T) %>%
        dplyr::summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))
    }else{
      data_stats <- data %>% dplyr::group_by(group_var, cluster, .drop=T) %>%
        dplyr::summarise(across(everything(), ~ median(.x, na.rm = TRUE)))
    }

    #data_stats <- data_stats %>% reshape2::melt(.,id.vars=c("cluster","group_var"))
    return(data_stats)
  }
}

#' PC_loadings
#'
#' @title PC_loadings
#' @description Function to visualize the effect each marker has on the principle component.
#' @param fcd flow cytometry dataset.
#' @param prefix Prefix for the output.
#' @param data_slot data slot to use for the calculation of the PCA, e.g. "orig" or "norm".
#' @param nPC Number of principle components to show.
#' @return Figure of PC loadings
#'
#' @import dplyr
#' @import cowplot
#' @import ggplot2
#'
#' @export


PC_loadings <- function(fcd,
                        prefix = NULL,
                        data_slot = "orig",
                        nPC = 3) {
  # select markers for PC_loading accordng to prefix (default -> all markers)
  markers <- used_markers(fcd = fcd, data_slot = data_slot, input_type = "pca", prefix = prefix, mute = T)

  pca_result <- prcomp(fcd$expr[[data_slot]][, colnames(fcd$expr[[data_slot]]) %in% markers, drop = F])
  pca_rotations <- as.matrix(pca_result$rotation)
  plot.list <- list()
  for (i in colnames(pca_rotations)) {
    tmp <- as.data.frame(pca_rotations)
    tmp$marker <- rownames(pca_rotations)
    tmp <- tmp[, c(i, "marker")]
    colnames(tmp) <- c("loadings", "marker")
    tmp <- tmp %>% dplyr::arrange(loadings)
    tmp$marker <- factor(tmp$marker, levels = tmp$marker)
    plot.list[[i]] <- ggplot(tmp, aes(x = loadings, y = marker)) +
      geom_point() +
      theme_linedraw() +
      ggtitle(paste(i,prefix))
  }
  cowplot::plot_grid(plotlist = plot.list[1:nPC], ncol = 3)
}


#' scaleColors
#'
#' @title scaleColors
#' @description Defines the color coding for an heatmap.
#' @param data data to use.
#' @param maxvalue value at which the color is fully red / blue.
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
#' @references code for this function was obtained from the R package ArchR (https://github.com/GreenleafLab/ArchR/blob/master/R/HelperUtils.R)
#'
#' @export
confusionMatrix <- function (i = NULL, j = NULL)
{
  ui <- unique(i)
  uj <- unique(j)
  m <- Matrix::sparseMatrix(i = match(i, ui), j = match(j,uj),
                            x = rep(1, length(i)), dims = c(length(ui), length(uj)))
  rownames(m) <- ui
  colnames(m) <- uj
  m
}

#' plot_dim_red
#'
#' @title Dimensionality reduction dotplot
#' @description \code{plot_dim_red()} generates a dotplot of the coordinates of any dimensionality reduction performed on a \code{condor} object. The plot can be colored by any variable both numeric (e.g. expression) or categorical (e.g. clustering/metadata).
#' @param fcd flow cytometry data set, that has been subjected to dimensionality reduction with cyCONDOR.
#' @param expr_slot expr_slot from which to take marker expression values, default is "orig".
#' Corrected input data should be handled cautiously.
#' @param reduction_method string specifying which dimensionality reduction method to use ("umap", "tSNE", "diffmap", "pca").
#' @param reduction_slot string specifying reduction name in reduction_method to use for visualization, e.g. "pca_orig".
#' @param cluster_slot string specifying which clustering slot to use to find variable specified in cluster_var.
#' @param add_pseudotime Logical, if plot should be colored by pseudotime.
#' @param pseudotime_slot string specifying pseudotime name to use for visualization.
#' @param param parameter to visualize in the plot, this can be either a continuous variable or a categorical one, the function will react differently accordingly.
#' @param order logical if you want to order the dots in the plot, by expression for example. This can help to find small populations of positive cells. If set to FALSE, the plotting order of the cells is randomized.
#' @param title title of the plot.
#' @param limX limits of the x axes (e.g. c(-1, 7)).
#' @param limY limits of the y axes (e.g. c(-1, 7)).
#' @param dot_size size of the dots.
#' @param alpha transparency of the dots.
#' @param color_discrete colors for discrete parameters, must be provided as vector of the same length as the number of factors of `param`.
#' @param color_gradient colors for continuous parameters.
#' @param remove_guide logical, if you want to remove the guide.
#' @param facet_by_variable option to facet the plot by a variable, if FALSE the plot is not faceted, if TRUE the plot is faceted by the `param` variable. If any other variable is provided (e.g. "group") the plot will be faceted by this variable.
#' @param label_clusters logical: If clusters should be labeled with a text box.
#' @param label_size size of the labels.
#' @param label_color color of the labels.
#' @param raster TRUE or FALSE, if plot should be returned as raster image, this option lowers the quality of the plot but makes it easier to work with images with high number of cells.
#' @param seed seed is set for reproducibility.
#' @import ggplot2
#' @import RColorBrewer
#' @import devtools
#' @import ggpubr
#' @import ggsci
#' @import ggrastr
#' @import ggrepel
#' @return plot marker or list of markers
#'
#' @export
plot_dim_red <- function(fcd,
                         expr_slot = NULL,
                         reduction_method,
                         reduction_slot,
                         cluster_slot = NULL,
                         add_pseudotime = FALSE,
                         pseudotime_slot,
                         param,
                         order = FALSE,
                         title = "Dimensionality Reduction Plot",
                         limX = NULL,
                         limY = NULL,
                         dot_size = 0.1,
                         alpha = 0.2,
                         color_discrete = cluster_palette,
                         color_gradient = colors,
                         remove_guide = FALSE,
                         facet_by_variable = FALSE,
                         label_clusters = FALSE,
                         label_size = 3.5,
                         label_color = "black",
                         raster = FALSE,
                         seed= 91) {

  # Preparing the data
  checkInput(fcd = fcd,
             check_expr_slot = F,
             check_cluster_slot = F,
             check_reduction = T,
             check_cell_anno = TRUE,
             expr_slot = expr_slot,
             cluster_slot = cluster_slot,
             reduction_method = reduction_method,
             reduction_slot = reduction_slot)

  data <- cbind(fcd[[reduction_method]][[reduction_slot]], fcd$anno$cell_anno)

  if (is.null(expr_slot) == FALSE) {

    checkInput(fcd = fcd,
               check_expr_slot = T,
               check_cluster_slot = F,
               check_reduction = T,
               expr_slot = expr_slot,
               cluster_slot = cluster_slot,
               reduction_method = reduction_method,
               reduction_slot = reduction_slot)

    data <- cbind(data, fcd$expr[[expr_slot]])

  }

  if (is.null(cluster_slot) == FALSE) {

    cluster_var <- "Description" # Bug fix, can be done better

    checkInput(fcd = fcd,
               check_expr_slot = F,
               check_cluster_slot = T,
               check_reduction = T,
               expr_slot = expr_slot,
               cluster_slot = cluster_slot,
               cluster_var = cluster_var,
               reduction_method = reduction_method,
               reduction_slot = reduction_slot)

    data <- cbind(data, fcd$clustering[[cluster_slot]])

  }

  if (isTRUE(add_pseudotime)) {

    data <- cbind(data, fcd$pseudotime[[pseudotime_slot]])

  }

  # Check if `param` exist
  if(!param %in% colnames(data)){
    stop('column "',param,'" is not available in the provided data')
  }

  # Check if `facet_by_variable` exist
  if (facet_by_variable != FALSE & facet_by_variable != TRUE) {

    if(!facet_by_variable %in% colnames(data)){
      stop('column "',facet_by_variable,'" is not available in the provided data')
    }

  }


  # Selection of raster plot or standard
  if (raster == FALSE) {

    core_numeric <- list(geom_point(aes(color = poi), alpha = alpha, size = dot_size),
                         theme_bw(),
                         theme(aspect.ratio = 1, panel.grid = element_blank()),
                         ggtitle(title),
                         scale_color_gradientn(colours = color_gradient),
                         labs(color = param))

    core_discrete <- list(geom_point(aes(color = poi), alpha = alpha, size = dot_size),
                          theme_bw(),
                          theme(aspect.ratio = 1, panel.grid = element_blank()),
                          ggtitle(title),
                          scale_colour_manual(values = color_discrete),
                          guides(color = guide_legend(override.aes = list(size=5, alpha = 1))),
                          labs(color = param))

  } else {

    core_numeric <- list(geom_point_rast(aes(color = poi), alpha = alpha, size = dot_size),
                         theme_bw(),
                         theme(aspect.ratio = 1, panel.grid = element_blank()),
                         ggtitle(title),
                         scale_color_gradientn(colours = color_gradient),
                         labs(color = param))

    core_discrete <- list(geom_point_rast(aes(color = poi), alpha = alpha, size = dot_size),
                          theme_bw(),
                          theme(aspect.ratio = 1, panel.grid = element_blank()),
                          ggtitle(title),
                          scale_colour_manual(values = color_discrete),
                          guides(color = guide_legend(override.aes = list(size=5, alpha = 1))),
                          labs(color = param))

  }


  if(is.numeric(data[[param]]) == TRUE){

    if(order == TRUE){
      data <- data[order(data[[param]], decreasing = F), ]
    }else if(order == FALSE){
      # order rows randomly for plotting
      set.seed(seed= seed)
      cells <- sample(x = rownames(data))
      data <- data[cells,]
    }

    colnames(data)[colnames(data) == param] <- "poi"

    # This if statement if to fix the faceting by selected variable issue with R4
    if (facet_by_variable != TRUE & facet_by_variable != FALSE) {

      colnames(data)[colnames(data) == facet_by_variable] <- "facet"

    }

    if (reduction_method == "umap") {

      p1 <- ggplot(data, aes(x = UMAP1, y = UMAP2)) + core_numeric

    }

    if (reduction_method == "diffmap") {

      p1 <- ggplot(data, aes(x = DC_1, y = DC_2)) + core_numeric

    }

    if (reduction_method == "tSNE") {

      p1 <- ggplot(data, aes(x = tSNE1, y = tSNE2)) + core_numeric

    }

    if (reduction_method == "pca") {

      p1 <- ggplot(data, aes(x = PC1, y = PC2)) + core_numeric

    }

    if (facet_by_variable != FALSE) {

      p1 <- p1 + facet_wrap(~facet)

    }

  }else {
    if(order == FALSE){
      # order rows randomly for plotting
      set.seed(seed= seed)
      cells <- sample(x = rownames(data))
      data <- data[cells,]
    }
    colnames(data)[colnames(data) == param] <- "poi"

    # This if statement is to fix the faceting by selected variable issue with R4
    if (facet_by_variable != TRUE & facet_by_variable != FALSE) {

      colnames(data)[colnames(data) == facet_by_variable] <- "facet"

    }

    if (reduction_method == "umap") {

      p1 <- ggplot(data, aes(x = UMAP1, y = UMAP2)) + core_discrete

    }

    if (reduction_method == "diffmap") {

      p1 <- ggplot(data, aes(x = DC_1, y = DC_2)) + core_discrete

    }

    if (reduction_method == "tSNE") {

      p1 <- ggplot(data, aes(x = tSNE1, y = tSNE2)) + core_discrete

    }

    if (reduction_method == "pca") {

      p1 <- ggplot(data, aes(x = PC1, y = PC2)) + core_discrete

    }

    if (facet_by_variable != FALSE) {

      if(facet_by_variable == TRUE) {

        p1 <- p1 + facet_wrap(~poi)

      } else {

        p1 <- p1 + facet_wrap(~facet)

      }

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

  if (label_clusters == TRUE) {

    data$poi_lab <- data$poi

    data$poi_lab <- ifelse(duplicated(data$poi_lab) == TRUE, NA, as.character(data$poi_lab))

    p1 <- p1 + geom_label_repel(aes(label = data$poi_lab), size = label_size, color = label_color, force = 10)

  }


  return(p1)

}


#' plot_marker_HM
#'
#' @title Heatmap of scaled expression by cell population
#' @description
#' \code{plot_marker_HM()} generates a heatmap of scaled mean marker expression for each cell population.
#' @param fcd flow cytometry data set, that has been subjected to clustering or cell type label prediction with cyCONDOR
#' @param expr_slot expr_slot from which to take marker expression values, default is "orig".
#' Corrected input data should be handled cautiously.
#' @param marker_to_exclude (optional) vector of characters indicating which features in expression matrix should not be included in the heatmap.
#' @param cluster_slot string specifying which clustering slot to use to find variable specified in cluster_var
#' @param cluster_var string specifying variable name in cluster_slot that identifies cell population labels to be used (e.g. clusters, metaclusters or predicted labels).
#' @param cluster_rows logical indicating if heatmap rows should be clustered (default: FALSE).
#' @param cluster_cols logical indicating if heatmap columns should be clustered (default: FALSE).
#' @param maxvalue max value for the coloring (default: NULL, automatically defined).
#' @param size size of the individual squares and font.
#' @param title character string, title of the plot
#' @import reshape2
#' @import Rmisc
#' @import pheatmap
#' @import Hmisc
#'
#' @returns
#' A heatmap of scaled mean expression, depicting markers in rows and cell populations in columns.
#'
#' @export
plot_marker_HM <- function(fcd,
                           expr_slot = "orig",
                           marker_to_exclude = NULL,
                           cluster_slot,
                           cluster_var,
                           cluster_rows = FALSE,
                           cluster_cols = FALSE,
                           maxvalue = NULL,
                           size = 10,
                           title = "Heatmap of scaled expression"
){

  #### check slots, cellIDs und varibles
  checkInput(fcd = fcd,
             check_expr_slot = T,
             check_cluster_slot = T,
             expr_slot = expr_slot,
             cluster_slot = cluster_slot,
             cluster_var = cluster_var)


  #### prepare data
  data <- fcd$expr[[expr_slot]]
  data$cluster <- fcd$clustering[[cluster_slot]][[cluster_var]]

  #### check if markers are present in expr slot
  if(!is.null(marker_to_exclude)){

    marker <- unique(marker_to_exclude)
    marker_present <- marker[marker %in% colnames(fcd$expr[[expr_slot]])]

    if(length(marker_present) < length(marker)){
      warning('The following markers could not be found in expr slot ',expr_slot,': ',
              paste(marker[!marker %in% marker_present], collapse = ","))
    }

    ## remove markers to exclude from data
    if(length(marker_present) > 0){
      data <- data[,!colnames(data) %in% marker_present]
    }

  }


  ## calculate mean
  tmp <- reshape2::melt(data, id.vars = c("cluster"))
  tmp <- Rmisc::summarySE(data = tmp, measurevar = "value", groupvars = c("cluster", "variable"))

  tmp <- reshape2::dcast(tmp[, c("cluster","variable","value")], variable ~ cluster)
  rownames(tmp) <- tmp$variable
  tmp$variable <- NULL

  ## scale values
  tmp <- t(base::scale(t(tmp)))

  #### plot
  p <- pheatmap::pheatmap(tmp, scale = "none",
                          breaks = cyCONDOR::scaleColors(data = tmp, maxvalue = maxvalue)[["breaks"]],
                          color = cyCONDOR::scaleColors(data = tmp, maxvalue = maxvalue)[["color"]],
                          main = title, cellwidth = size, cellheight = size,
                          cluster_rows = cluster_rows, cluster_cols = cluster_cols)

  return(p)
}


#' plot_confusion_HM
#'
#' @title plot_confusion_HM
#' @description \code{plot_confusion_HM()} generates a heatmap showing the contribution of each group_var to a cell population after normalizing all levels in group_var to the same cell numbers.
#' @param fcd flow cytometry data set, that has been subjected to clustering or cell type label prediction with cyCONDOR
#' @param cluster_slot string specifying which clustering slot to use to find variable specified in cluster_var
#' @param cluster_var string specifying variable name in cluster_slot that identifies cell population labels to be used to calculate the confusion.
#' @param group_var string indicating variable name in cell_anno that should be used to calculate the relative contribution to the variable specified in cluster_var.
#' @param numeric logical, indicating if levels in cluster_var should be ordered in ascending numerical order.
#' @param size size of the individual squares and font
#' @param title character string, title of the plot.
#' @param cluster_cols logical indicating if columns should be clustered (default: FALSE)
#' @param cluster_rows logical indicating if rows should be clustered (default: FALSE)
#' @return \code{plot_confusion_HM()} first calculates cell counts for each combination of group_var and cell population and normalizes the counts to a total of 1000 cells per group_var.
#' Afterwards the percentage of cells coming from each level in group_var is calculated per cell population. The normalization of counts corrects the visualization for differences in total cells (events) measured per group_var.
#' @import pheatmap
#'
#' @export
plot_confusion_HM <- function(fcd,
                              cluster_slot,
                              cluster_var,
                              group_var,
                              numeric = FALSE,
                              size = 15,
                              title = "Confusion matrix",
                              cluster_cols = FALSE,
                              cluster_rows = FALSE) {

  #### check slots, cellIDs und varibles
  checkInput(fcd = fcd,
             check_cluster_slot = T,
             check_cell_anno = T,
             cluster_slot = cluster_slot,
             cluster_var = cluster_var,
             group_var = group_var)


  #### prepare data
  data <- data.frame(cellID = rownames(fcd$clustering[[cluster_slot]]),
                     cluster = fcd$clustering[[cluster_slot]][[cluster_var]],
                     group_var = fcd$anno$cell_anno[[group_var]])


  ## cell counts per cluster-group_var combination
  cells_cluster <- cyCONDOR::confusionMatrix(paste0(data$cluster),
                                             paste0(data$group_var))


  if (numeric == TRUE) {
    cells_cluster <- cells_cluster[order(as.numeric(rownames(cells_cluster))),order(colnames(cells_cluster))]
  }else{
    cells_cluster <- cells_cluster[,order(colnames(cells_cluster))]
  }

  cells_cluster <- as.matrix(cells_cluster)

  ## normalize to 1000 cells per group_var
  tmp <- round(t(t(cells_cluster)/colSums(cells_cluster))*1000,3)
  ## calculate percentage of cells from group_var per cluster
  scaled_cM <- round((tmp / Matrix::rowSums(tmp))*100,2)

  p <- pheatmap::pheatmap(
    mat = t(scaled_cM),
    border_color = "black",display_numbers = TRUE,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    cellwidth = size,
    cellheight = size,
    main = title)

  return(p)
}

#' plot_frequency_boxplot
#'
#' @title box plot of cell population frequencies
#' @description \code{plot_frequency_boxplot()} plots cell population frequencies of samples grouped by group_var as boxplots.
#' @param fcd flow cytometry data set, that has been subjected to clustering or cell type label prediction with cyCONDOR
#' @param cluster_slot string specifying which clustering slot to use to find variable specified in cluster_var
#' @param cluster_var string specifying variable name in cluster_slot that identifies cell population labels to be used (e.g. clusters, metaclusters or predicted labels).
#' @param groups_to_show vector of strings indicating levels in group_var that should be included for plotting. By default all groups are plotted.
#' @param group_var string indicating variable name in cell_anno that should be used to group samples in sample_var
#' @param sample_var string indicating variable name in cell_anno that defines sample IDs to be used.
#' @param numeric logical, if TRUE numeric levels in cluster_var are ordered in increasing order and "Cluster_" is pasted before number, if FALSE alphabetical ordering is applied.
#' @param color_palette vector of colors to be used to fill box plots
#' @param dot_size Size of the dots.
#' @import dplyr
#' @import reshape2
#' @import ggplot2
#' @returns \code{plot_frequency_boxplots()} returns a list of boxplots, each element of the list contains the plot of one cell population.
#'
#' @export
plot_frequency_boxplot<-function(fcd,
                                 cluster_slot,
                                 cluster_var,
                                 sample_var,
                                 group_var,
                                 groups_to_show = NULL,
                                 numeric = F,
                                 color_palette = cluster_palette,
                                 dot_size = 2){

  #### check slots, cell IDs and variables
  checkInput(fcd = fcd,
             check_cluster_slot = T,
             check_cell_anno = T,
             cluster_slot = cluster_slot,
             cluster_var = cluster_var,
             group_var = group_var,
             sample_var = sample_var)


  #### prepare data
  data<-data.frame(cellID=rownames(fcd$clustering[[cluster_slot]]),
                   cluster=fcd$clustering[[cluster_slot]][[cluster_var]],
                   group_var=fcd$anno$cell_anno[[group_var]],
                   sample_var=fcd$anno$cell_anno[[sample_var]])

  if(is.numeric(data$group_var)){
    data$group_var <- as.character(data$group_var)}
  if(is.numeric(data$sample_var)){
    data$sample_var <- as.character(data$sample_var)}

  ## prepare frequency table
  tmp <- cyCONDOR::confusionMatrix(paste0(data$sample_var), paste0(data$cluster))
  tmp <- as.matrix(tmp)
  if (numeric == TRUE) {
    tmp <- tmp[order(rownames(tmp)), order(as.numeric(colnames(tmp)))]
    colnames(tmp) <- paste0("Cluster_", colnames(tmp))
  } else if (numeric == FALSE) {
    tmp <- tmp[order(rownames(tmp)), order(colnames(tmp))]
  } else {
    #tmp <- tmp[order(rownames(tmp)), order(colnames(tmp))]
    stop('argument "numeric" needs to be set to TRUE or FALSE.')
  }

  tmp <- tmp/rowSums(tmp) * 100
  tmp <- as.data.frame(tmp)
  tmp$sample_var <- rownames(tmp)
  tmp <- dplyr::left_join(tmp, unique(data[,c("group_var","sample_var")]), by = "sample_var")

  ## check if samples are uniquely assigned to one group
  anno <- unique(data[,c("sample_var","group_var")])
  if(!nrow(anno) == length(unique(anno$sample_var))){
    stop("levels in sample_var cannot be uniquely asigned to one group_var. Make sure that each sample_var is unambiguous associated with only one group in group_var.")
  }


  ## select groups to be included in plotting
  if (!is.null(groups_to_show)) {
    tmp <- base::subset(tmp, group_var %in% groups_to_show)

    if(nrow(tmp)==0){
      stop('none of groups spedified in "groups_to_show" are present in argument "group_var"')
    }
  }


  #### plot
  tmp<-reshape2::melt(tmp, id.vars = c("group_var", "sample_var"))

  plot.list <- list()
  for (i in unique(tmp$variable)) {
    data_cluster <- tmp[tmp$variable ==i,]

    p <- ggplot(data_cluster, aes(x = group_var, y = value, fill = group_var)) +
      geom_boxplot(outlier.colour = NA) +
      geom_point(shape = 21, fill = "white", color = "black", size = dot_size, position = position_dodge(0.4)) +
      theme_bw() +
      theme(#aspect.ratio = 2,
        panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle=90,hjust = 1, vjust = 0.5)) +
      ggtitle(i) +
      expand_limits(y = 0) + xlab("") + ylab("percentage") +
      scale_fill_manual(values = color_palette)

    plot.list[[i]] <- p
    rm(p, data_cluster)
  }

  return(plot.list)
}

#' plot_frequency_barplot
#'
#' @title bar chart of cell population frequencies
#' @description \code{plot_frequency_barplot()} plots cell population frequencies for each level in group_var as stacked barplot.
#' @param fcd flow cytometry data set, that has been subjected to clustering or cell type label prediction with cyCONDOR
#' @param cluster_slot string specifying which clustering slot to use to find variable specified in cluster_var.
#' @param cluster_var string specifying variable name in cluster_slot that identifies cell population labels to be used (e.g. clusters, metaclusters or predicted labels).
#' @param group_var string indicating variable name in cell_anno that defines grouping variable to be used (x-axis), e.g. group or sample ID.
#' @param facet_var (optional) string indicating variable name in cell_anno that should be used to group levels in group_var via faceting.
#' @param color_palette vector of colors to be used to fill bar chart plots
#' @param title title of the plot, default is "Frequency"
#' @import ggplot2
#' @returns \code{plot_frequency_barplot} returns a plot showing cell population frequencies for each level in group_var. The plot is faceted by another variable, when provided a facet_var.
#'
#' @export
plot_frequency_barplot<-function (fcd = condor,
                                  cluster_slot,
                                  cluster_var,
                                  group_var,
                                  facet_var = NULL,
                                  color_palette = cluster_palette,
                                  title = "Frequency") {

  #### check slots, cell IDs and variables
  checkInput(fcd = fcd,
             check_cluster_slot = T,
             check_cell_anno = T,
             cluster_slot = cluster_slot,
             cluster_var = cluster_var,
             group_var = group_var)

  if(!is.null(facet_var)){
    if(!facet_var %in% colnames(fcd$anno$cell_anno)){
      stop('column "',facet_var,'" is not available in cell_anno')
    }
  }


  #### prepare data
  data <- data.frame(cellID=rownames(fcd$clustering[[cluster_slot]]),
                     cluster=fcd$clustering[[cluster_slot]][[cluster_var]],
                     group_var=fcd$anno$cell_anno[[group_var]])

  if(!is.null(facet_var)){
    data$facet_var <- fcd$anno$cell_anno[[facet_var]]
  }

  if(is.numeric(data$cluster)){
    data$cluster <- as.character(data$cluster)
  }
  if(is.numeric(data$group_var)){
    data$group_var <- as.character(data$group_var)
  }

  #### plot
  p <- ggplot(data, aes(x = group_var, fill = cluster)) +
    geom_bar(position = "fill", color = "black", linewidth = 0.3) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    labs(y = "Frequency (%)", x = "", fill = cluster_var) + ggtitle(title) +
    scale_fill_manual(values = color_palette)

  if(!is.null(facet_var)==T){
    p <- p + facet_wrap(.~facet_var, scales = "free_x")
  }

  return(p)
}

#' plot_counts_barplot
#'
#' @title bar chart of cell population counts
#' @description \code{plot_counts_barplot()} plots cell population counts for each level in group_var as bar plot.
#' @param fcd flow cytometry data set, that has been subjected to clustering or cell type label prediction with cyCONDOR
#' @param cluster_slot string specifying which clustering slot to use to find variable specified in cluster_var
#' @param cluster_var string specifying variable name in cluster_slot that identifies cell population labels to be used (e.g. clusters, metaclusters or predicted labels).
#' @param group_var string indicating variable name in cell_anno that defines grouping variable to be used (x-axis), e.g. group or sample ID.
#' @param facet_var (optional) string indicating variable name in cell_anno that should be used to group levels in group_var via faceting.
#' @param facet_by_clustering logical, if set to TRUE, the plot will be faceted by cell populations specified in cluster_var. If set to FALSE not faceting by population will be applied (default). If both, facet_var and facet_by_clustering are used plot will be faceted by both using applying facet_grid.
#' @param facet_ncol numeric, indicating how many columns faceted plot should have, if either facet_var or Facet_by_clustering is used.
#' @param color_palette vector of colors to be used to fill bar chart plots
#' @param title title of the plot, default is "Counts"
#' @import ggplot2
#' @returns \code{plot_counts_barplot()} returns a plot showing cell numbers per cell population for each level in group_var. The plot is faceted by another variable, when provided a facet_var and or by the cell population, depending on facet_by_clustering.
#' @export
plot_counts_barplot<-function (fcd = condor,
                               cluster_slot,
                               cluster_var,
                               group_var,
                               facet_var = NULL,
                               facet_by_clustering = F,
                               facet_ncol = NULL,
                               color_palette = cluster_palette,
                               title = "Counts") {

  #### check slots, cell IDs and variables
  checkInput(fcd = fcd,
             check_cluster_slot = T,
             check_cell_anno = T,
             cluster_slot = cluster_slot,
             cluster_var = cluster_var,
             group_var = group_var)

  if(!is.null(facet_var)){
    if(!facet_var %in% colnames(fcd$anno$cell_anno)){
      stop('column "',facet_var,'" is not available in cell_anno')
    }
  }
  if(!(isTRUE(facet_by_clustering) || isFALSE(facet_by_clustering))){
    stop('argument "facet_by_clustering" needs to be TRUE or FALSE.')
  }


  #### prepare data
  data <- data.frame(cellID=rownames(fcd$clustering[[cluster_slot]]),
                     cluster=fcd$clustering[[cluster_slot]][[cluster_var]],
                     group_var=fcd$anno$cell_anno[[group_var]])

  if(!is.null(facet_var)){
    data$facet_var<-fcd$anno$cell_anno[[facet_var]]
  }

  if(is.numeric(data$cluster)){
    data$cluster <- as.character(data$cluster)
  }
  if(is.numeric(data$group_var)){
    data$group_var <- as.character(data$group_var)
  }

  #### plot
  p <- ggplot(data, aes(x = group_var, fill = cluster)) +
    geom_bar(position = "stack", color = "black", linewidth = 0.3) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    labs(y = "cell count", x = "", fill = cluster_var) + ggtitle(title) +
    scale_fill_manual(values = color_palette)

  if(!is.null(facet_var) & isFALSE(facet_by_clustering)){
    p <- p + facet_wrap(.~facet_var, scales = "free_x", ncol = facet_ncol)
  }
  if(is.null(facet_var) & isTRUE(facet_by_clustering)){
    p <- p + facet_wrap(.~cluster, scales = "free_x", ncol = facet_ncol)+
      theme(legend.position = "none")
  }
  if(!is.null(facet_var) & isTRUE(facet_by_clustering)){
    p <- p + facet_grid(cluster~facet_var, scales = "free")+
      theme(legend.position = "none")
  }

  return(p)
}

#' plot_dim_density
#'
#' @title dimensionality reduction plot with density distribution
#' @description \code{plot_dim_density()} plots the density distribution of each level of a grouping variable on top of dimensionality reduction plot. At the moment, only UMAP visualization is supported
#' @param fcd flow cytometry data set comprising dimensionality reduction calculated with cyCONDOR.
#' @param reduction_method string specifying which dimensionality reduction method to use. At the moment, only "umap" is supported.
#' @param reduction_slot string specifying reduction name in reduction_method to use for visualization, e.g. "pca_orig".
#' @param group_var string indicating variable name in cell_anno for which density distribution will be plotted for each level in group_var.
#' @param title character string, title of the plot
#' @param dot_size size of the background dots.
#' @param alpha transparency of the background dots.
#' @param color_density vector of strings indicating the color palette from \code{RColorBrewer} to be used for the density gradient. A names vector using levels in group_var allows the assignment of a specific color to a group level. If vector is not names, entries will be assigned to groups in alphabetic order and color palettes will be reused if vector has less colors than levels in group_var
#' @import ggrastr
#' @import ggpubr
#' @import ggplot2
#'
#' @return Density Plot

#'
#' @export
plot_dim_density <- function(fcd,
                             reduction_method = "umap",
                             reduction_slot,
                             group_var,
                             title = "",
                             dot_size = 0.1,
                             alpha = 0.2,
                             color_density = c("Blues", "Oranges", "Reds")) {



  #### check slots, cellIDs und varibles
  checkInput(fcd = fcd,
             check_cell_anno = T,
             check_reduction = T,
             reduction_method = reduction_method,
             reduction_slot = reduction_slot,
             group_var = group_var)


  #### prepare data
  data <- as.data.frame(fcd[[reduction_method]][[reduction_slot]])
  data$group_var <- fcd$anno$cell_anno[[group_var]]

  groups <- unique(data$group_var)
  groups <- groups[order(groups)]

  #### assign colors, if vector is not named or adjusted
  if(is.null(names(color_density))){

    if(length(color_density) < length(groups)){
      color_density<-rep_len(color_density,length(groups))
    }
    names(color_density) <- groups
  }


  #### plotting
  if (reduction_method == "umap"){

    p0 <- ggplot(data = data, aes(x = UMAP1, y = UMAP2)) +
      ggrastr::geom_point_rast(fill="grey",color="grey",size=dot_size,alpha=alpha) +
      theme_bw() +
      theme(aspect.ratio = 1, panel.grid = element_blank())

    plot.list <- list()
    for (i in groups) {

      p <- p0 +
        stat_density_2d(data = data[data$group_var == i,],
                        #aes(x = UMAP1, y = UMAP2, fill = ..level..),
                        aes(x = UMAP1, y = UMAP2, fill = after_stat(level)),
                        geom = "polygon", contour = T, h = .5) +
        scale_fill_distiller(palette = color_density[i], direction = 1) +
        ggtitle(i)

      plot.list[[i]] <- p
    }
  }


  out <- ggpubr::ggarrange(plotlist = plot.list, legend = "right", ncol = 3, common.legend = FALSE)
  out <- ggpubr::annotate_figure(out, top = ggpubr::text_grob(title))

  return(out)

}

#' plot_marker_group_HM
#'
#' @title Heatmap of scaled expression to compare two groups
#' @description
#' \code{plot_marker_group_HM()} generates a heatmap of scaled mean marker expression for each cell population split by a grouping variable.
#' @param fcd flow cytometry data set, that has been subjected to clustering or cell type label prediction with cyCONDOR
#' @param expr_slot expr_slot from which to take marker expression values, default is "orig".
#' Corrected input data should be handled cautiously.
#' @param marker_to_exclude (optional) vector of characters indicating which features in expression matrix should not be included in the heatmap.
#' @param cluster_slot string specifying which clustering slot to use to find variable specified in cluster_var
#' @param cluster_var string specifying variable name in cluster_slot that identifies cell population labels to be used (e.g. clusters, metaclusters or predicted labels).
#' @param group_var string indicating variable name in cell_anno that should be used for subgrouping each cell population.
#' @param maxvalue max value for the coloring (default: NULL, automatically defined).
#' @param size size of the individual squares and font.
#' @param title character string, title of the plot
#' @import reshape2
#' @import Rmisc
#' @import pheatmap
#'
#' @returns
#' A heatmap of scaled mean expression, depicting markers in rows and cell populations grouped by a grouping variable in columns.
#'
#' @export
plot_marker_group_HM <- function(fcd,
                                 expr_slot = "orig",
                                 marker_to_exclude = NULL,
                                 cluster_slot,
                                 cluster_var,
                                 group_var,
                                 maxvalue = NULL,
                                 size = 10,
                                 title = "Heatmap of scaled expression"
){

  #### check slots, cellIDs und varibles
  checkInput(fcd = fcd,
             check_expr_slot = T,
             check_cluster_slot = T,
             check_cell_anno = T,
             expr_slot = expr_slot,
             cluster_slot = cluster_slot,
             cluster_var = cluster_var,
             group_var = group_var)


  #### prepare data
  data <- fcd$expr[[expr_slot]]
  data$cluster <- fcd$clustering[[cluster_slot]][[cluster_var]]
  data$group_var <- fcd$anno$cell_anno[[group_var]]


  #### check if markers are present in expr slot
  if(!is.null(marker_to_exclude)){

    marker <- unique(marker_to_exclude)
    marker_present <- marker[marker %in% colnames(fcd$expr[[expr_slot]])]

    if(length(marker_present) < length(marker)){
      warning('The following markers could not be found in expr slot ',expr_slot,': ',
              paste(marker[!marker %in% marker_present], collapse = ","))
    }

    ## remove markers to exclude from data
    if(length(marker_present) > 0){
      data <- data[,!colnames(data) %in% marker_present]
    }
  }


  ## calculate mean
  tmp <- reshape2::melt(data, id.vars = c("cluster","group_var"))
  tmp <- Rmisc::summarySE(data = tmp, measurevar = "value", groupvars = c("cluster", "variable", "group_var"))
  ngroups <- length(unique(tmp$group_var))

  tmp <- reshape2::dcast(tmp[, c("cluster", "variable", "group_var","value")], variable ~ cluster + group_var)
  rownames(tmp) <- tmp$variable
  tmp$variable <- NULL

  ## scale values
  tmp <- t(base::scale(t(tmp)))

  #### plot
  p <- pheatmap::pheatmap(tmp, scale = "none",
                          cluster_rows = FALSE, cluster_cols = FALSE,
                          breaks = cyCONDOR::scaleColors(data = tmp,maxvalue = maxvalue)[["breaks"]],
                          color = cyCONDOR::scaleColors(data = tmp, maxvalue = maxvalue)[["color"]],
                          main = title,
                          cellwidth = size, cellheight = size, fontsize = size,
                          gaps_col = seq(0, ncol(tmp), ngroups))

  return(p)
}

#' plot_marker_violinplot
#'
#' @title violin plot of marker expression in cell populations
#' @description
#' \code{plot_marker_violinplot()} plots the expression of selected markers as violin plots for cell populations.
#' @param fcd flow cytometry data set, that has been subjected to clustering or cell type label prediction with cyCONDOR
#' @param marker vector of characters indicating which features of the expression matrix should be plotted.
#' @param expr_slot expr_slot from which to take marker expression values, default is "orig".
#' Corrected input data should be handled cautiously.
#' @param cluster_slot string specifying which clustering slot to use to find variable specified in cluster_var
#' @param cluster_var string specifying variable name in cluster_slot that identifies cell population labels to be used (e.g. clusters, metaclusters or predicted labels).
#' @param cluster_to_show vector of strings indicating levels in cluster_var that should be included for plotting.
#' @param group_var (optional) string indicating variable name in cell_anno that should be used to split violin plots.
#' @param color_palette vector of colors to be used to fill violin plots, when group_var is used
#' @import ggplot2
#' @returns \code{plot_marker_violinplot()} returns either one plot in case only one marker is provided via \code{marker} argument or a list of plots, if several markers are requested.
#' @details The violin plots are plotted with default parameters of ggplot2's \code{\link[ggplot2]{geom_violin}} and horizontal lines indicate the median.
#'
#' @export
#'
plot_marker_violinplot<- function(fcd,
                                  marker,
                                  expr_slot = "orig",
                                  cluster_slot,
                                  cluster_var,
                                  cluster_to_show = NULL,
                                  group_var = NULL,
                                  color_palette = cluster_palette){

  #### check slots, cellIDs und varibles
  checkInput(fcd = fcd,
             check_expr_slot = T,
             check_cluster_slot = T,
             check_cell_anno = T,
             expr_slot = expr_slot,
             cluster_slot = cluster_slot,
             cluster_var = cluster_var,
             group_var = group_var)


  #### check if markers are present in expr slot
  marker <- unique(marker)
  marker_present <- marker[marker %in% colnames(fcd$expr[[expr_slot]])]
  if(length(marker_present) == 0){
    stop('None of the requested markers is present in expr slot "',expr_slot,'".')
  }
  if(length(marker_present) < length(marker)){
    warning('The following marker could not be found in expr slot ',expr_slot,': ',paste(marker[!marker %in% marker_present], collapse = ", "))
  }


  #### get cluster of interest
  if(!is.null(cluster_to_show)){
    ## remove potential duplicates
    cluster_to_show <- unique(cluster_to_show)

    ## get cluster of interest
    cluster_present <- cluster_to_show[cluster_to_show %in% unique(fcd$clustering[[cluster_slot]][[cluster_var]])]
    if(length(cluster_present) == 0){
      stop('None of the provided clusters are present in selected cluster_var "',cluster_var,'".')
    }
    if(length(cluster_present) < length(cluster_to_show)){
      warning('The following clusters could not be found in cluster_var "',cluster_var,'": ',
              paste(cluster_to_show[!cluster_to_show %in% cluster_present], collapse = ","))
    }
  }else{
    ## default
    cluster_present <- unique(fcd$clustering[[cluster_slot]][[cluster_var]])
  }


  #### prepare data
  data <- fcd$expr[[expr_slot]]
  data$cluster <- fcd$clustering[[cluster_slot]][[cluster_var]]


  if(is.numeric(data$cluster)){
    data$cluster <- as.character(data$cluster)
  }
  if(!is.null(group_var)){
    data$group_var <- fcd$anno$cell_anno[[group_var]]

    if(is.numeric(data$group_var)){
      data$group_var <- as.character(data$group_var)
    }
  }

  ## subset to cluster of interest
  data <- data[data$cluster %in% cluster_present,]


  #### plotting
  plot.list<-list()

  if(!is.null(group_var)){

    for (i in marker_present){
      data_sub<-data[,c(as.character(i),"cluster","group_var")]
      colnames(data_sub)[colnames(data_sub) == as.character(i)] <- "marker"

      p <- ggplot(data_sub, aes(cluster, marker,fill = group_var))+
        geom_violin(draw_quantiles = 0.5, trim = T, scale = "area")+
        theme_linedraw()+
        ggtitle(i)+
        theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5))+
        labs(x = cluster_var, y = "expression", fill = group_var) +
        scale_fill_manual(values = color_palette)

      plot.list[[i]] <- p
    }

  }else{

    for (i in marker_present){
      data_sub <- data[,c(as.character(i),"cluster")]
      colnames(data_sub)[colnames(data_sub) == as.character(i)] <- "marker"

      p <- ggplot(data_sub, aes(cluster, marker))+
        geom_violin(draw_quantiles = 0.5, trim = T, scale = "area", fill = "grey80")+
        theme_linedraw()+
        ggtitle(i)+
        theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5))+
        labs(x = cluster_var, y = "expression")

      plot.list[[i]] <- p
    }

  }

  #### return
  if(length(plot.list) == 1){
    p <- plot.list[[1]]
    return(p)

  }else{
    # p<-cowplot::plot_grid(plotlist = plot.list, ncol = 1)
    # return(p)
    return(plot.list)
  }

}


#' plot_marker_boxplot
#'
#' @title plot box plots of marker median or mean
#' @description \code{plot_marker_boxplot()} generates faceted box plots of aggregated expression per sample grouped by a grouping of interest for each marker-cluster combination.
#' @param fcd flow cytometry data set, that has been subjected to clustering or cell type label prediction with cyCONDOR
#' @param marker marker in provided expr_slot that should be included in plotting. By default all features are plotted.
#' @param expr_slot expr_slot from which to take marker expression values, default is "orig".
#' Corrected input data should be handled cautiously.
#' @param cluster_slot string specifying which clustering slot to use to find variable specified in cluster_var
#' @param cluster_var string specifying variable name in cluster_slot that identifies cell population labels to be used (e.g. clusters, metaclusters or predicted labels).
#' @param cluster_to_show vector of characters indicating levels in cluster_var that should be included for plotting.
#' By default all groups are plotted.
#' @param group_var string indicating variable name in cell_anno that should be used to group samples in sample_var
#' @param sample_var string indicating variable name in cell_anno that defines sample IDs to be used.
#' @param fun string indicating which aggregation of expression should be performed for each sample. Needs to be one out of "median" (default) or "mean".
#' @param facet_by_clustering logical indicating how to structure faceting of the plot. If set to FALSE (default) plot will be faceted by markers and clusters (or cell populations) are shown on x-axis. If set to TRUE, it is the other way around.
#' @param facet_ncol Number of columns to be used for faceting.
#' @param color_palette vector of colors that should be used to fill box plots.
#' @param dot_size numeric indicating the size of the dots.
#' @import dplyr
#' @import reshape2
#' @import ggplot2
#' @returns
#' A faceted ggplot showing box plots of aggregated expression per sample grouped by group_var.
#'
#' @export
plot_marker_boxplot<- function(fcd,
                               marker = NULL,
                               expr_slot ="orig",
                               cluster_slot = "Phenograph_pca_orig_k_60",
                               cluster_var = "metaclusters",
                               cluster_to_show = NULL,
                               group_var,
                               sample_var,
                               fun = "median",
                               facet_by_clustering = F,
                               facet_ncol = 5,
                               color_palette = cluster_palette,
                               dot_size=2){

  #### check slots, cellIDs und varibles
  checkInput(fcd = fcd,
             check_expr_slot = T,
             check_cluster_slot = T,
             check_cell_anno = T,
             expr_slot = expr_slot,
             cluster_slot = cluster_slot,
             cluster_var = cluster_var,
             group_var = group_var,
             sample_var = sample_var)

  if(!(isTRUE(facet_by_clustering) || isFALSE(facet_by_clustering))){
    stop('argument "facet_by_clustering" needs to be TRUE or FALSE.')
  }


  #### get markers of interest
  if(!is.null(marker)){
    # remove potential duplicates
    marker <- unique(marker)

    ## check if markers are present in expr slot
    marker_present <- marker[marker %in% colnames(fcd$expr[[expr_slot]])]
    if(length(marker_present) == 0){
      stop('None of the provided markers are present in expr slot "',expr_slot,'".')
    }
    if(length(marker_present) < length(marker)){
      warning('The following markers could not be found in expr slot "',expr_slot,'": ',
              paste(marker[!marker %in% marker_present], collapse = ","))
    }
  }else{
    ## default
    marker_present <- colnames(fcd$expr[[expr_slot]])
  }

  #### get cluster of interest
  if(!is.null(cluster_to_show)){
    ## remove potential duplicates
    cluster_to_show <- unique(cluster_to_show)

    ## get cluster of interest
    cluster_present <- cluster_to_show[cluster_to_show %in% unique(fcd$clustering[[cluster_slot]][[cluster_var]])]
    if(length(cluster_present) == 0){
      stop('None of the provided clusters are present in selected cluster_var "',cluster_var,'".')
    }
    if(length(cluster_present) < length(cluster_to_show)){
      warning('The following clusters could not be found in cluster_var "',cluster_var,'": ',
              paste(cluster_to_show[!cluster_to_show %in% cluster_present], collapse = ","))
    }
  }else{
    ## default
    cluster_present <- unique(fcd$clustering[[cluster_slot]][[cluster_var]])
  }


  #### prepare data
  data <- fcd$expr[[expr_slot]]
  data$cluster <- fcd$clustering[[cluster_slot]][[cluster_var]]
  data$group_var <- fcd$anno$cell_anno[[group_var]]
  data$sample_var <- fcd$anno$cell_anno[[sample_var]]

  if(is.numeric(data$cluster)){
    data$cluster <- as.character(data$cluster)
  }
  if(is.numeric(data$group_var)){
    data$group_var <- as.character(data$group_var)
  }

  ## subset to marker of interest
  data <- as.data.frame(data[,c(marker_present,"cluster","group_var","sample_var")])

  ## subset to cluster of interest
  data <- data[data$cluster %in% cluster_present,]

  ## check if samples are uniquely assigned to one group
  anno <- unique(data[,c("sample_var","group_var")])
  if(!nrow(anno) == length(unique(anno$sample_var))){
    stop("levels in sample_var cannot be uniquely asigned to one group_var. Make sure that each sample_var is unambiguous associated with only one group in group_var.")
  }

  #### calculate median or mean per sample for each marker-cluster combination
  if(fun == "mean"){
    data_stats <- data %>% dplyr::group_by(sample_var, group_var, cluster, .drop=T) %>%
      dplyr::summarise(dplyr::across(dplyr::everything(), ~ mean(.x, na.rm = TRUE)))
  }else if(fun == "median"){
    data_stats <- data %>% dplyr::group_by(sample_var, group_var, cluster, .drop=T) %>%
      dplyr::summarise(dplyr::across(dplyr::everything(), ~ median(.x, na.rm = TRUE)))
  }else{
    stop('aggregation function "fun" must be either set to "mean" or to "median".')
  }

  data_stats <- data_stats %>% reshape2::melt(.,id.vars=c("cluster","group_var","sample_var"))
  n_groups <- length(unique(data_stats$group_va))
  #### plot
  if(facet_by_clustering == F){
    p <- ggplot(data_stats, aes(x = cluster, y = value, fill = group_var)) +
      geom_boxplot(outlier.colour = NA,position = position_dodge2(0.85, preserve = "single")) +
      geom_point(position=position_jitterdodge(jitter.width = 0.05,jitter.height = 0),
                 shape = 21, color = "black", size = dot_size) +
      theme_bw() +
      theme(#aspect.ratio = 2,
        panel.grid = element_blank(),
        legend.position = "right",
        axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) +
      expand_limits(y = 0) + xlab(cluster_var) + ylab(paste(fun," expression")) + labs(fill = group_var) +
      facet_wrap(.~variable, scales = "free_y", ncol = facet_ncol) +
      scale_fill_manual(values = color_palette)

  }else if(facet_by_clustering == T){
    p <- ggplot(data_stats, aes(x = variable, y = value, fill = group_var)) +
      geom_boxplot(outlier.colour = NA, position = position_dodge2(0.85, preserve = "single")) +
      geom_point(position = position_jitterdodge(jitter.width = 0.05, jitter.height = 0),
                 shape = 21, color = "black", size = dot_size) +
      theme_bw() +
      theme(#aspect.ratio = 2,
        panel.grid = element_blank(),
        legend.position = "right",
        axis.text.x = element_text(angle=90,hjust = 1, vjust = 0.5)) +
      expand_limits(y = 0) + xlab("marker") + ylab(paste(fun," expression")) + labs(fill = group_var) +
      facet_wrap(.~cluster, scales = "free_y", ncol = facet_ncol)+
      scale_fill_manual(values = color_palette)

  }else{
    stop('argument "facet_by_clustering" needs to be a logical.')
  }

  return(p)
}

#' plot_marker_dotplot
#'
#' @title plot_marker_dotplot
#' @description
#' \code{plot_marker_dotplot()} generates a classical scatter plot of two markers
#' @param fcd flow cytometry data set, that has been subjected to clustering or cell type label prediction with cyCONDOR
#' @param expr_slot expr_slot from which to take marker expression values, default is "orig".
#' Corrected input data should be handled cautiously.
#' @param marker_x marker name in expr slot of which transformed expression is shown on x-axis.
#' @param marker_y marker name in expr slot of which transformed expression is shown on y-axis.
#' @param cluster_slot string specifying which clustering slot to use to find variable specified in cluster_var.
#' @param cluster_var string specifying variable name in cluster_slot that identifies cell population labels to be used (e.g. clusters, metaclusters or predicted labels).
#' @param cluster_to_show (optional) vector of characters indicating levels in cluster_var that should be included for plotting.
#' By default all groups are plotted.
#' @param group_var (optional) string indicating variable name in cell_anno that should be used to color dots in plot. If not used, dots are clustered by cluster_var.
#' @param order logical if you want to order the dots in the plot by the variable used for coloring (TRUE). If set to FALSE (default), order of the cells is randomized.
#' @param seed a seed is set for reproducibility of the plotting result.
#' @param color_palette vector of colors that should be used to color dots.
#' @param dot_size numeric indicating the size of the dots.
#' @param title Title of the plot.
#' @import ggplot2
#' @return The function returns a scatter plot of two features available in the expression matrix. By default, dots are colored by cell population label provided in cluster_var. If coloring by a metadata is wanted instead, a group_var can be defined. Further, if only a selection of levels available in cluster_var should be included in plotting, a vector of labels of interest can be provided to cluster_to_show argument.
#'
#'
#' @export
plot_marker_dotplot <- function(fcd,
                                expr_slot = "orig",
                                marker_x,
                                marker_y,
                                cluster_slot,
                                cluster_var,
                                cluster_to_show = NULL,
                                group_var = NULL,
                                order = F,
                                seed = 91,
                                color_palette = cluster_palette,
                                dot_size=2,
                                title ="") {

  #### check slots, cellIDs und varibles
  checkInput(fcd = fcd,
             check_expr_slot = T,
             check_cluster_slot = T,
             check_cell_anno = T,
             expr_slot = expr_slot,
             cluster_slot = cluster_slot,
             cluster_var = cluster_var,
             group_var = group_var)


  #### check markers
  if(!marker_x %in% colnames(fcd$expr[[expr_slot]])){
    stop('marker "',marker_x, '" is not present in expr slot "',expr_slot,'".')
  }
  if(!marker_y %in% colnames(fcd$expr[[expr_slot]])){
    stop('marker "',marker_y, '" is not present in expr slot "',expr_slot,'".')
  }


  #### get cluster of interest
  if(!is.null(cluster_to_show)){
    ## remove potential duplicates
    cluster_to_show <- unique(cluster_to_show)

    ## get cluster of interest
    cluster_present <- cluster_to_show[cluster_to_show %in% unique(fcd$clustering[[cluster_slot]][[cluster_var]])]
    if(length(cluster_present) == 0){
      stop('None of the provided clusters are present in selected cluster_var "',cluster_var,'".')
    }
    if(length(cluster_present) < length(cluster_to_show)){
      warning('The following clusters could not be found in cluster_var "',cluster_var,'": ',
              paste(cluster_to_show[!cluster_to_show %in% cluster_present], collapse = ","))
    }
  }else{
    ## default
    cluster_present <- unique(fcd$clustering[[cluster_slot]][[cluster_var]])
  }


  #### prepare data
  data <- fcd$expr[[expr_slot]]
  data$cluster <- fcd$clustering[[cluster_slot]][[cluster_var]]

  if(is.numeric(data$cluster)){
    data$cluster <- as.character(data$cluster)
  }
  if(!is.null(group_var)){
    data$group_var <- fcd$anno$cell_anno[[group_var]]

    if(is.numeric(data$group_var)){
      data$group_var <- as.character(data$group_var)
    }
  }

  ## select marker of interest
  colnames(data)[colnames(data) == marker_x] <- "X"
  colnames(data)[colnames(data) == marker_y] <- "Y"

  ## subset to cluster of interest
  data <- data[data$cluster %in% cluster_present,]


  #### plot

  if(!is.null(group_var)){

    ## order of dots
    if(order == T){
      #data <- data[order(data$group_var, decreasing = F), ]
    }else if(order == F){
      # order rows randomly for plotting
      set.seed(seed = seed)
      cells <- sample(x = rownames(data))
      data <- data[cells,]
    }else{
      stop('argument order needs to be FALSE or TRUE.')
    }

    p <- ggplot(data, aes(x = X, y = Y)) +
      geom_point(aes(color = group_var), size = dot_size) +
      theme_bw()+
      theme(aspect.ratio = 1, panel.grid = element_blank())+
      ggtitle(title) +
      scale_colour_manual(values = color_palette) +
      guides(color = guide_legend(override.aes = list(size=5, alpha = 1))) +
      expand_limits(y = 0) + xlab(marker_x) + ylab(marker_y) + labs(fill = group_var)
  }else{

    ## order of dots
    if(order == T){
      #data <- data[order(data$cluster, decreasing = F), ]
    }else if(order == F){
      # order rows randomly for plotting
      set.seed(seed= seed)
      cells <- sample(x = rownames(data))
      data <- data[cells,]
    }else{
      stop('argument order needs to be FALSE or TRUE.')
    }

    p <- ggplot(data, aes(x = X, y = Y)) +
      geom_point(aes(color = cluster), size = dot_size) +
      theme_bw()+
      theme(aspect.ratio = 1, panel.grid = element_blank())+
      ggtitle(title) +
      scale_colour_manual(values = color_palette) +
      guides(color = guide_legend(override.aes = list(size=5, alpha = 1))) +
      expand_limits(y = 0) + xlab(marker_x) + ylab(marker_y) + labs(fill = cluster_var)
  }

  return(p)
}

#' plot_marker_density
#'
#' @title density plot of marker expression in cell populations grouped by a meta variable
#' @description
#' \code{plot_marker_density()} plots the density distribution of expression for cell populations. Each line indicates distribution for one level in group_var.
#' @param fcd flow cytometry data set, that has been subjected to clustering or cell type label prediction with cyCONDOR
#' @param marker vector of characters indicating which features of the expression matrix should be plotted.
#' @param expr_slot expr_slot from which to take marker expression values, default is "orig".
#' Corrected input data should be handled cautiously.
#' @param cluster_slot string specifying which clustering slot to use to find variable specified in cluster_var
#' @param cluster_var string specifying variable name in cluster_slot that identifies cell population labels to be used (e.g. clusters, metaclusters or predicted labels).
#' @param cluster_to_show (optional) vector of strings indicating levels in cluster_var that should be included for plotting.
#' @param group_var string indicating variable name in cell_anno for which density distributions get calculated.
#' @param facet_var (optional) string indicating variable name in cell_anno that should be used to facet the plots by a meta variable.
#' @param facet_ncol numeric, number of columns used for faceting.
#' @param color_palette vector of colors to be used to color the density lines
#' @returns \code{plot_marker_density()} returns either one plot in case only one marker is provided via \code{marker} argument or a list of plots, if several markers are requested.
#' @details The density plots are plotted with default parameters of ggplot2's \code{\link[ggplot2]{geom_density}}.
#' @import ggplot2
#'
#' @export
plot_marker_density <- function(fcd,
                                marker,
                                expr_slot = "orig",
                                cluster_slot,
                                cluster_var,
                                cluster_to_show = NULL,
                                group_var,
                                facet_var = NULL,
                                facet_ncol = 5,
                                color_palette = cluster_palette) {

  #### check slots, cellIDs und varibles
  checkInput(fcd = fcd,
             check_expr_slot = T,
             check_cluster_slot = T,
             check_cell_anno = T,
             expr_slot = expr_slot,
             cluster_slot = cluster_slot,
             cluster_var = cluster_var,
             group_var = group_var)


  #### check if markers are present in expr slot
  marker <- unique(marker)
  marker_present <- marker[marker %in% colnames(fcd$expr[[expr_slot]])]
  if(length(marker_present) == 0){
    stop('None of the requested markers is present in expr slot "',expr_slot,'".')
  }
  if(length(marker_present) < length(marker)){
    warning('The following marker could not be found in expr slot ',expr_slot,': ',paste(marker[!marker %in% marker_present], collapse = ", "))
  }


  #### get cluster of interest
  if(!is.null(cluster_to_show)){

    ## remove potential duplicates
    cluster_to_show <- unique(cluster_to_show)

    ## get cluster of interest
    cluster_present <- cluster_to_show[cluster_to_show %in% unique(fcd$clustering[[cluster_slot]][[cluster_var]])]
    if(length(cluster_present) == 0){
      stop('None of the provided clusters are present in selected cluster_var "',cluster_var,'".')
    }
    if(length(cluster_present) < length(cluster_to_show)){
      warning('The following clusters could not be found in cluster_var "',cluster_var,'": ',
              paste(cluster_to_show[!cluster_to_show %in% cluster_present], collapse = ","))
    }
  }else{
    ## default
    cluster_present <- unique(fcd$clustering[[cluster_slot]][[cluster_var]])
  }


  #### prepare data
  data <- fcd$expr[[expr_slot]]
  data$cluster <- fcd$clustering[[cluster_slot]][[cluster_var]]
  data$group_var <- fcd$anno$cell_anno[[group_var]]

  ## add optional facet_var
  if (!is.null(facet_var)) {
    if(facet_var %in% colnames(fcd$anno$cell_anno)){
      data$facet_var <- fcd$anno$cell_anno[[facet_var]]
    }else{
      stop('facet_var: column "',facet_var,'" is not available in cell_anno')
    }
  }
  ## fix numeric group_var
  if(is.numeric(data$group_var)){
    data$group_var <- as.character(data$group_var)
  }

  ## subset to cluster of interest
  data <- data[data$cluster %in% cluster_present,]


  #### plotting
  plot.list<-list()

  if(!is.null(facet_var)){
    for (i in marker_present){
      data_sub <- data[,c(as.character(i),"cluster","group_var","facet_var")]
      colnames(data_sub)[colnames(data_sub) == as.character(i)] <- "marker"

      p <- ggplot(data_sub, aes(x = marker)) +
        geom_density(aes(color = group_var)) +
        theme_bw() +
        scale_colour_manual(values = color_palette) +
        theme(aspect.ratio = 1, panel.grid = element_blank()) +
        labs(color = group_var) +
        ggtitle(i) + xlab("expression") +
        facet_grid(facet_var~cluster)

      plot.list[[i]] <- p
    }
  }else{
    for (i in marker_present){
      data_sub <- data[,c(as.character(i),"cluster","group_var")]
      colnames(data_sub)[colnames(data_sub) == as.character(i)] <- "marker"

      p <- ggplot(data_sub, aes(x = marker)) +
        geom_density(aes(color = group_var)) +
        theme_bw() +
        scale_colour_manual(values = color_palette) +
        theme(aspect.ratio = 1, panel.grid = element_blank()) +
        labs(color = group_var) +
        ggtitle(i) + xlab("expression") +
        facet_wrap(.~cluster, ncol = facet_ncol)

      plot.list[[i]] <- p
    }
  }

  #### return
  if(length(plot.list) == 1){
    p <- plot.list[[1]]
    return(p)

  }else{
    return(plot.list)
  }
}

#' plot_marker_ridgeplot
#'
#' @title ridgeline plot of marker expression in cell populations
#' @description
#' \code{plot_marker_ridgeplot()} plots the expression of selected markers as ridgeline plots for cell populations.
#' @param fcd flow cytometry data set, that has been subjected to clustering or cell type label prediction with cyCONDOR
#' @param marker vector of characters indicating which features of the expression matrix should be plotted.
#' @param expr_slot expr_slot from which to take marker expression values, default is "orig".
#' Corrected input data should be handled cautiously.
#' @param cluster_slot string specifying which clustering slot to use to find variable specified in cluster_var
#' @param cluster_var string specifying variable name in cluster_slot that identifies cell population labels to be used (e.g. clusters, metaclusters or predicted labels).
#' @param cluster_to_show (optional) vector of strings indicating levels in cluster_var that should be included for plotting.
#' @param color_palette vector of colors to be used to fill density distributions.
#' @param alpha numeric, adjust alpha to be used on fill color of density distributions.
#' @import ggridges
#' @returns \code{plot_marker_ridgeplot()} returns either one plot in case only one marker is provided via \code{marker} argument or a list of plots, if several markers are requested.
#' @details The ridgeline plots are plotted with default parameters of ggridges's \code{\link[ggridges]{geom_density_ridges}}
#' @export
#'
plot_marker_ridgeplot<- function(fcd,
                                 marker,
                                 expr_slot = "orig",
                                 cluster_slot,
                                 cluster_var,
                                 cluster_to_show = NULL,
                                 color_palette = cluster_palette,
                                 alpha = 1){

  #### check slots, cellIDs und varibles
  checkInput(fcd = fcd,
             check_expr_slot = T,
             check_cluster_slot = T,
             check_cell_anno = T,
             expr_slot = expr_slot,
             cluster_slot = cluster_slot,
             cluster_var = cluster_var)


  #### check if markers are present in expr slot
  marker <- unique(marker)
  marker_present <- marker[marker %in% colnames(fcd$expr[[expr_slot]])]
  if(length(marker_present) == 0){
    stop('None of the requested markers is present in expr slot "',expr_slot,'".')
  }
  if(length(marker_present) < length(marker)){
    warning('The following marker could not be found in expr slot ',expr_slot,': ',paste(marker[!marker %in% marker_present], collapse = ", "))
  }


  #### get cluster of interest
  if(!is.null(cluster_to_show)){
    ## remove potential duplicates
    cluster_to_show <- unique(cluster_to_show)

    ## get cluster of interest
    cluster_present <- cluster_to_show[cluster_to_show %in% unique(fcd$clustering[[cluster_slot]][[cluster_var]])]
    if(length(cluster_present) == 0){
      stop('None of the provided clusters are present in selected cluster_var "',cluster_var,'".')
    }
    if(length(cluster_present) < length(cluster_to_show)){
      warning('The following clusters could not be found in cluster_var "',cluster_var,'": ',
              paste(cluster_to_show[!cluster_to_show %in% cluster_present], collapse = ","))
    }
  }else{
    ## default
    cluster_present <- unique(fcd$clustering[[cluster_slot]][[cluster_var]])
  }


  #### prepare data
  data <- fcd$expr[[expr_slot]]
  data$cluster <- fcd$clustering[[cluster_slot]][[cluster_var]]

  ## subset to cluster of interest
  data <- data[data$cluster %in% cluster_present,]

  if(is.numeric(data$cluster)){
    data$cluster <- as.character(data$cluster)
  }


  #### plotting
  plot.list<-list()

  for (i in marker_present){
    data_sub <- data[,c(as.character(i),"cluster")]
    colnames(data_sub)[colnames(data_sub) == as.character(i)] <- "marker"

    if(is.factor(data_sub$cluster)){
      data_sub$cluster <- factor(data_sub$cluster, levels = rev(levels(data_sub$cluster)))
    }

    p <- ggplot(data_sub, aes(x = marker, y = cluster, fill = cluster)) +
      ggridges::geom_density_ridges(alpha = alpha) +
      ggridges::theme_ridges(center_axis_labels = T) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_discrete(expand = expand_scale(mult = c(0.01, 0.2))) +
      scale_fill_manual(values = color_palette) +
      theme(legend.position = "none") +
      xlab("expression") + ylab(cluster_var) + ggtitle(i)

    plot.list[[i]] <- p
  }


  #### return
  if(length(plot.list) == 1){
    p <- plot.list[[1]]
    return(p)

  }else{
    return(plot.list)
  }

}
