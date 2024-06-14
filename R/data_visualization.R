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
  markers <- used_markers(fcd = fcd, data_slot = data_slot, method = "pca", prefix = prefix, mute = T)

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

#' plot_marker
#'
#' @title plot_marker
#' @description Function to visualize a selected parameter on a dimensionality reduction.
#' @param data This is the input dataframe for the visualization, this part of the code can still be improved but for the moment you have to cbind the dataframe with the informations you want to plot.
#' @param param Parameter to visualize in the plot, this can be either a continous variable or a categorical one, the function will react differently according.
#' @param order Logical if you want to order the dots in the plot, by expression for example. This can help to find small populations of positive cells. If set to FALSE, the plotting order of the cells is randomized.
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
#' @param label_clusters Logical: If clusters should be labeled with a text box.
#' @param label_size Size of the labels.
#' @param label_color Color of the labels.
#' @param raster TRUE or FALSE, if plot should be returned as raster image
#' @param seed A seed is set for reproducibility.
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
                            facet_by_variable = FALSE,
                            label_clusters = FALSE,
                            label_size = 3.5,
                            label_color = "black",
                            raster = FALSE,
                            seed= 91) {

  # Selection of raster plot or standard
  if (raster == FALSE) {

    core_numeric <- list(geom_point(aes(color = poi), alpha = apha, size = dot_size),
                         theme_bw(),
                         theme(aspect.ratio = 1, panel.grid = element_blank()),
                         ggtitle(title),
                         scale_color_gradientn(colours = color_gradient),
                         labs(color = param))

    core_discrete <- list(geom_point(aes(color = poi), alpha = apha, size = dot_size),
                          theme_bw(),
                          theme(aspect.ratio = 1, panel.grid = element_blank()),
                          ggtitle(title),
                          scale_colour_manual(values = color_discrete),
                          guides(color = guide_legend(override.aes = list(size=5, alpha = 1))),
                          labs(color = param))

  } else {

    core_numeric <- list(geom_point_rast(aes(color = poi), alpha = apha, size = dot_size),
                         theme_bw(),
                         theme(aspect.ratio = 1, panel.grid = element_blank()),
                         ggtitle(title),
                         scale_color_gradientn(colours = color_gradient),
                         labs(color = param))

    core_discrete <- list(geom_point_rast(aes(color = poi), alpha = apha, size = dot_size),
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

    if (dim_red == "UMAP") {

      p1 <- ggplot(data, aes(x = UMAP1, y = UMAP2)) + core_numeric

    }

    if (dim_red == "DM") {

      p1 <- ggplot(data, aes(x = DC_1, y = DC_2)) + core_numeric

    }

    if (dim_red == "tSNE") {

      p1 <- ggplot(data, aes(x = tSNE1, y = tSNE2)) + core_numeric

    }

    if (dim_red == "PCA") {

      p1 <- ggplot(data, aes(x = PC1, y = PC2)) + core_numeric

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

    if (dim_red == "UMAP") {

      p1 <- ggplot(data, aes(x = UMAP1, y = UMAP2)) + core_discrete

    }

    if (dim_red == "DM") {

      p1 <- ggplot(data, aes(x = DC_1, y = DC_2)) + core_discrete

    }

    if (dim_red == "tSNE") {

      p1 <- ggplot(data, aes(x = tSNE1, y = tSNE2)) + core_discrete

    }

    if (dim_red == "PCA") {

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
#' `plot_marker_HM()` generates a heatmap of scaled mean marker expression for each cell population.
#' @param fcd flow cytometry data set, that has been subjected to clustering or cell type label prediction with *cyCONDOR*
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
#' @description `plot_confusion_HM()` generates a heatmap showing the contribution of each group_var to a cell population after normalizing all levels in group_var to the same cell numbers.
#' @param fcd flow cytometry data set, that has been subjected to clustering or cell type label prediction with *cyCONDOR*
#' @param cluster_slot string specifying which clustering slot to use to find variable specified in cluster_var
#' @param cluster_var string specifying variable name in cluster_slot that identifies cell population labels to be used to calculate the confusion.
#' @param group_var string indicating variable name in cell_anno that should be used to calculate the relative contribution to the variable specified in cluster_var.
#' @param numeric logical, indicating if levels in cluster_var should be ordered in ascending numerical order.
#' @param size size of the individual squares and font
#' @param title character string, title of the plot.
#' @param cluster_cols logical indicating if columns should be clustered (default: FALSE)
#' @param cluster_rows logical indicating if rows should be clustered (default: FALSE)
#' @return `plot_confusion_HM()` first calculates cell counts for each combination of group_var and cell population and normalizes the counts to a total of 1000 cells per group_var.
#' Afterwards the percentage of cells coming from each level in group_var is calculated per cell population. The normalization of counts corrects the visualization for differences in total cells (events) measured per group_var.
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

#' boxplot_and_stats
#'
#' @title boxplot_and_stats
#' @description Calculate boxplot and summary statistics for a wanted comparison.
#' @param annotation Sample annotation to be used for the plot.
#' @param sample_var Column name containing the sample IDs.
#' @param group_var Column name defining the groupping for plotting.
#' @param which_groups Grouping to be used.
#' @param variable Variable used to stratify the plotting.
#' @param numeric Logical if the groupping is numeric.
#' @param test.type Test to be performed. (see need some development here).
#' @param paired.test Logicle: If the test should be paired.
#' @return boxplot_and_stats
#'
#' @export
boxplot_and_stats <- function(annotation,
                              sample_var,
                              group_var,
                              which_groups = NULL,
                              variable,
                              numeric = TRUE,
                              test.type,
                              paired.test = FALSE) {

  container <- list()

  tmp <- confusionMatrix(paste0(annotation[[sample_var]]),
                         paste0(variable))

  tmp <- as.matrix(tmp)

  if (numeric == TRUE) {

    tmp <- tmp[order(rownames(tmp)), order(as.numeric(colnames(tmp)))]

  }

  if (numeric == FALSE) {

    tmp <- tmp[order(rownames(tmp)), order(colnames(tmp))]

  }

  colnames(tmp) <- paste("Cluster_", colnames(tmp), sep = "")

  tmp <- tmp/rowSums(tmp)*100

  tmp <- as.data.frame(tmp)

  tmp[[sample_var]] <- rownames(tmp)

  tmp <- merge(tmp, unique(annotation), by = sample_var)

  container[["per_data"]] <- tmp

  colnames(tmp)[colnames(tmp) == group_var] <- "groups"

  if (!is.null(which_groups)) {

    tmp <- subset(tmp, groups %in% which_groups)

  }

  for (variable in colnames(tmp)[grep("Cluster_", colnames(tmp))]) {

    data <- tmp

    colnames(data)[colnames(data) == variable] <- "poi"

    p1 <- ggplot(data, aes(x = groups, y = poi, fill = groups))+
      geom_boxplot(outlier.colour = NA)+
      geom_point(shape = 21, fill = "white", color = "black",
                 size = 4, position = position_dodge(0.4))+
      theme_bw()+ theme(aspect.ratio = 2, panel.grid = element_blank())+
      #scale_fill_manual(values = cluster_palette)+
      ggtitle(variable)+
      expand_limits(y = 0)+
      xlab("")+
      stat_compare_means(method = test.type,
                         paired = paired.test,
                         label = "p.signif",
                         label.x.npc = 0.5)

    container[["plot"]][[variable]] <- p1

    rm(p1, data)

  }

  tmp2 <- suppressMessages(melt(tmp[, c(colnames(tmp)[grep("Cluster_", colnames(tmp))], "groups")]))

  stat <- compare_means(formula = value~groups,
                        data = tmp2,
                        paired = paired.test,
                        method = test.type,
                        group.by = "variable",
                        p.adjust.method = "none")

  container[["stats"]] <- stat

  return(container)

}

#' barplot_frequency
#'
#' @title barplot_frequency
#' @description This function output a stacked barplot for the cellular frequencies.
#' @param x_axes Groupping of the x axes.
#' @param colour Stratification to use on the stacked barplot.
#' @param color_palette color_palette.
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

#' plot_density
#'
#' @title plot_density
#' @description Density Plot of dimensionality reduction.
#' @param data data to be used for plotting.
#' @param param Parameter to be used for the density splitting.
#' @param title Title of the plot.
#' @param dim_red Dimensionality reduction to be used.
#' @param dot_size Size of the background dots.
#' @param alpha Transparency of the background dots.
#' @param color_density Colors of the density maps.
#' @return Density Plot
#'
#' @export
plot_density <- function(data,
                         param,
                         title = "adjust the title",
                         dim_red,
                         dot_size = 0.1,
                         alpha = 0.2,
                         color_density = c("Blues", "Oranges", "Reds")) {

  group_density <- list()

  colnames(data)[colnames(data) == param] <- "poi"

  sequence <- unique(data$poi)

  sequence <- sequence[order(sequence)]

  names(color_density) <- sequence

  if (dim_red == "UMAP"){
    p0 <- ggplot(data = data, aes(x = UMAP1, y = UMAP2))
  }

  if (dim_red == "DM"){
    p0 <- ggplot(data = data, aes(x = DM_1, y = DM_2))
  }

  for (i in sequence) {

    p <- p0 +
      geom_point_rast(fill="grey",color="grey",size=dot_size,alpha=alpha) +
      theme_bw() +
      theme(aspect.ratio = 1, panel.grid = element_blank()) +
      stat_density_2d(data = data[data$poi == i,],
                      aes(x = UMAP1, y = UMAP2, fill = ..level..),
                      geom = "polygon", contour = T, h = .5) +
      scale_fill_distiller(palette = color_density[i], direction = 1) +
      ggtitle(paste(dim_red, i, sep = "_"))

    group_density[[i]] <- p

  }

  out <- ggarrange(plotlist = group_density, legend = "right", ncol = 3, common.legend = FALSE)

  out <- annotate_figure(out, top = text_grob(title))

  return(out)

}

#' plot_marker_group_HM
#'
#' @title Heatmap of scaled expression to compare two groups
#' @description
#' `plot_marker_group_HM()` generates a heatmap of scaled mean marker expression for each cell population split by a grouping variable.
#' @param fcd flow cytometry data set, that has been subjected to clustering or cell type label prediction with *cyCONDOR*
#' @param expr_slot expr_slot from which to take marker expression values, default is "orig".
#' Corrected input data should be handled cautiously.
#' @param marker_to_exclude (optional) vector of characters indicating which features in expression matrix should not be included in the heatmap.
#' @param cluster_slot string specifying which clustering slot to use to find variable specified in cluster_var
#' @param cluster_var string specifying variable name in cluster_slot that identifies cell population labels to be used (e.g. clusters, metaclusters or predicted labels).
#' @param group_var string indicating variable name in cell_anno that should be used for subgrouping each cell population.
#' @param maxvalue max value for the coloring (default: NULL, automatically defined).
#' @param size size of the individual squares and font.
#' @param title character string, title of the plot
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

#' violinplot_marker
#'
#' @title violinplot_marker
#' @description Calculate a violin plot for the expression of selected markers.
#' @param fcd flow cytometry dataset.
#' @param data_slot data to use for the calculation, e.g. "expr" or "pca".
#' @param cluster_method methods for the clustering to use.
#' @param cluster_type type of clustering (slot name).
#' @param marker marker to plot.
#' @import Hmisc
#' @return violinplot_marker
#'
#' @export
violinplot_marker <- function(fcd,
                              data_slot = "orig",
                              cluster_method,
                              cluster_type,
                              marker = colnames(df[,-which(names(df) %in% c("metaclusters", "Description",
                                                                            "Phenograph", "FlowSOM"))])) {

  df <- cbind(fcd$clustering[[cluster_method]], fcd$expr[[data_slot]])

  plot.list <- list()

  for (i in intersect(colnames(df), marker)) {

    df.short <- df[, c(cluster_type, i)]
    colnames(df.short) <- c("cluster_type", "marker")

    plot.list[[i]] <- ggplot(df.short, aes(x = cluster_type, y = marker))+
      geom_violin()+
      xlab("Cluster")+
      ylab(i)+
      theme_linedraw()+
      stat_summary(fun.data = mean_sdl, size = 0.2,
                   geom = "pointrange", color = "red")+
      ggtitle(i)
  }

  cowplot::plot_grid(plotlist = plot.list, ncol = 1)
}

#' plot_marker_dotplot
#'
#' @title plot_marker_dotplot
#' @description
#' `plot_marker_dotplot` generates a classical scatter plot of two markers
#' @param fcd flow cytometry data set, that has been subjected to clustering or cell type label prediction with *cyCONDOR*
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
#' @return
#' The function returns a scatter plot of two features available in the expression matrix. By default, dots are colored by cell population label provided in cluster_var. If coloring by a metadata is wanted instead, a group_var can be defined. Further, if only a selection of levels available in cluster_var should be included in plotting, a vector of labels of interest can be provided to cluster_to_show argument.
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

  if(!is.null(group_var)){
    data$group_var <- fcd$anno$cell_anno[[group_var]]
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
      data <- data[order(data$group_var, decreasing = F), ]
    }else if(order == F){
      # order rows randomly for plotting
      set.seed(seed= seed)
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
      data <- data[order(data$cluster, decreasing = F), ]
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

#' densityplot_marker
#'
#' @title dotplot_cyto
#' @description Plot the distribution of the expression of selected markers, similar to FlowJo Histogram.
#' @param data Data to use for the plot.
#' @param marker Marker to be visualized.
#' @param split_by Variable used to split the plot.
#' @param color_by Variable used to color the plot.
#' @param color_discrete Color palette (default: cluster_palette)
#' @param title Plot title.
#' @return densityplot marker.
#'
#' @export
densityplot_marker <- function(data,
                               marker,
                               split_by = NULL,
                               color_by = NULL,
                               color_discrete = cluster_palette,
                               title) {

  colnames(data)[colnames(data) == marker] <- "poi"

  core <- list(geom_density(),
               theme_bw(),
               theme(aspect.ratio = 1, panel.grid = element_blank()),
               ggtitle(title))

  if (!is.null(color_by)) {

    colnames(data)[colnames(data) == color_by] <- "color"

    core <- list(geom_density(aes(color = color)),
                 theme_bw(),
                 scale_colour_manual(values = color_discrete),
                 theme(aspect.ratio = 1, panel.grid = element_blank()),
                 labs(color = color_by),
                 ggtitle(title),
                 xlab(marker))

  }

  if (!is.null(split_by)) {

    colnames(data)[colnames(data) == split_by] <- "split"
  }

  p1 <- ggplot(data, aes(x = poi)) + core

  if (!is.null(split_by)) {

    p1 <- p1 + facet_wrap(~split)

  }

  return(p1)

}
