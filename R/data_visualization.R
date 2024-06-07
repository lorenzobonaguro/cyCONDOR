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
#' @param seed A seed is set for reproducibility of the plotting result.
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
plot_marker_new <- function(data,
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
                            seed= 42) {

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
#' @param cluster_rows Logical: If heatmap rows should be clustered.
#' @param cluster_cols Logical: If heatmap columns should be clustered.
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
                     cluster_rows = cluster_rows, cluster_cols = cluster_cols)

}

#' confusion_HM
#'
#' @title confusion_HM
#' @description Visualize the confusion heatmap.
#' @param variables Variable to by used to calculate the confusion.
#' @param group Groupping to calculate the relative contribution to the variable.
#' @param size Size of the heatmap blocks.
#' @param title Title of the plot.
#' @param cluster_cols TRUE or FALSE if columns should be clustered
#' @param cluster_rows TRUE or FALSE if rows should be clustered
#' @return confusion_HM
#'
#' @export
confusion_HM <- function(variables,
                         group,
                         size = 15,
                         title = "Choose a good title for the plot",
                         cluster_cols = FALSE,
                         cluster_rows = FALSE) {

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
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
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

#' HM_differential_marker
#'
#' @title HM_differential_marker
#' @description Calculate HM to compare marker expression in different groups.
#' @param fcd flow cytometry dataset.
#' @param data_slot data to use for the calculation of the UMAP, e.g. "expr" or "pca".
#' @param cluster_method Name of the clustering result to be used.
#' @param cluster_type Type of the clustering.
#' @param maxvalue Max value for the coloring (Default: NULL, automatically defined).
#' @param size Size of the individual squares.
#' @param title Title for the plot.
#' @param exclusion Vector with markers to exclude from the plot.
#' @param group_by Grouping variable for visualization.
#' @return df_frequency
#'
#' @export
HM_differential_marker <- function(fcd,
                                   data_slot = "orig",
                                   cluster_method, #"Phenograph_pca_norm_k60", "FlowSOM_pca_orig_k5" or others
                                   cluster_type, #"Phenograph", "FlowSOM", "metaclusters"
                                   maxvalue = NULL,
                                   group_by, #"group"
                                   size = 10,
                                   title = "Choose a good title for the plot",
                                   exclusion = c("SSC-A", "FSC-A")) {

  tmp <- suppressMessages(melt(cbind(fcd$expr[[data_slot]],
                                     fcd$clustering[[cluster_method]],
                                     fcd$anno$cell_anno[colnames(fcd$anno$cell_anno) == group_by])))

  colnames(tmp)[colnames(tmp) == group_by] <- "group"

  #summarySE function from Rmisc package
  tmp <- summarySE(data = tmp, measurevar = "value", groupvars = c(cluster_type, "variable", "group"))

  colnames(tmp)[colnames(tmp) == cluster_type] <- "cluster"

  ngroups <- length(unique(tmp$group))

  tmp <- dcast(tmp[, c(1,2,3,5)], variable ~ cluster + group)

  rownames(tmp) <- tmp$variable #change rownames to marker names
  tmp$variable <- NULL #remove column "variable" with marker names
  tmp <- tmp[!rownames(tmp) %in% exclusion,] #exclude marker from df/plot

  tmp <- t(scale(t(tmp)))

  pheatmap::pheatmap(tmp, scale = "none", cluster_rows = FALSE, cluster_cols = FALSE,
                     breaks = scaleColors(data = tmp, maxvalue = maxvalue)[["breaks"]],
                     color = scaleColors(data = tmp, maxvalue = maxvalue)[["color"]],
                     main = title, cellwidth = size, cellheight = size,
                     gaps_col = seq(0, ncol(tmp), ngroups)) #creates gap between clusters, not sure if needed but might be helpful for a big heatmap

}

#' PC_loading
#'
#' @title PC_loading
#' @description Calculate PCA loading.
#' @param fcd flow cytometry dataset.
#' @param data_slot data to use for the calculation of the UMAP, e.g. "expr" or "pca".
#' @param number Number of principal components to plot.
#' @return PC loadings
#'
#' @export
PC_loadings <- function(fcd, data_slot = "orig", number) {

  pca_result <- prcomp(fcd[["expr"]][[data_slot]])

  pca_rotations <- as.matrix(pca_result$rotation)

  plot.list <- list()

  for (i in colnames(pca_rotations)){

    tmp <- as.data.frame(pca_rotations)
    tmp$marker <- rownames(pca_rotations)
    tmp <- tmp[,c(i,"marker")]
    colnames(tmp) <- c("loadings", "marker")
    tmp <- tmp %>% dplyr::arrange(loadings)
    tmp$marker <- factor(tmp$marker, levels=tmp$marker)

    plot.list[[i]] <- ggplot(tmp,aes(x=loadings,y=marker))+
      geom_point()+
      theme_linedraw()+
      ggtitle(i)
  }

  cowplot::plot_grid(plotlist = plot.list[1:number], ncol = 3)

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

#' dotplot_cyto
#'
#' @title dotplot_cyto
#' @description Plot a classical dotplot such as normally seen in FC analysis (e.g. FlowJo).
#' @param data data to use for the plot.
#' @param subset LOGICLE, if TRUE the data will be subsetted for selected population(s).
#' @param annotation Vector of the same lenght at data containing the variable to subset.
#' @param color_by Parameter to be used to define the plot colors.
#' @param subset_char Vector of variables to keep in the plot (e.g. c("CD4T, "CD8T")).
#' @param color_discrete Color palette to use (default: cluster palette).
#' @param x Parameter to display in the X axes.
#' @param y Parameter to display in the Y axes.
#' @param title Plot title.
#' @return 2D dotplot.
#'
#' @export
dotplot_cyto <- function(data,
                         subset = FALSE,
                         annotation,
                         color_by,
                         subset_char,
                         color_discrete = cluster_palette,
                         x,
                         y,
                         title) {

  if (subset == TRUE) {

    data <- data[annotation %in% subset_char,]

    color_by <- color_by[annotation %in% subset_char]

  }

  colnames(data)[colnames(data) == x] <- "X"

  colnames(data)[colnames(data) == y] <- "Y"

  p1 <- ggplot(data, aes(x = X, y = Y)) +
    geom_point(aes(color = color_by)) +
    theme_bw()+
    theme(aspect.ratio = 1, panel.grid = element_blank())+
    ggtitle(title) +
    scale_colour_manual(values = color_discrete) +
    guides(color = guide_legend(override.aes = list(size=5, alpha = 1))) +
    xlab(x) +
    ylab(y)

  return(p1)

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
