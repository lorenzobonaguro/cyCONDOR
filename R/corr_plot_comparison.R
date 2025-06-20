# This function is designed to calculate a sample-wise correlation of cell population percentages that are found in CyCondor and FlowJo.


#' corr_plot_comparison
#' @title Sample-wise correlation of cell type proportions to compare manual gating and cyCONDOR.
#' @description 'corr_plot_comparison' performs a sample-wise correlation of cell type proportions obtained via manual gating and via clustering on cyCONDOR. The result is shown in a correlogram.
#' @param condor_df data frame containing the cell type frequencies obtained via clustering and annotation in cyCONDOR.
#' @param flowjo_df data frame containing the cell type frequencies obtained via manual gating and annotation (e.g. FlowJo).
#' @param sample_col name of the column containing sample names. This column name needs to be matching between the two data frames.
#' @param method_corr correlation method used in the correlogram (default = "pearson").



corr_plot_comparison <- function(condor_df,
                                 flowjo_df,
                                 sample_col,
                                 method_corr = "pearson"){


     # Merge the two dataframes by Sample/DONOR ID
     merged_data <- merge(condor_df, flowjo_df, by = sample_col)

     flowjo_cols <- colnames(merged_data)[grepl("FlowJo",colnames(merged_data))]
     condor_cols <- colnames(merged_data)[grepl("Condor",colnames(merged_data))]


     condor_clean <- setNames(condor_cols, sub("^[^_]+_", "", condor_cols))
     flowjo_clean <- setNames(flowjo_cols, sub("^[^_]+_", "", flowjo_cols))


     matching_cell_types <- intersect(names(condor_clean), names(flowjo_clean))

     total_cell_types_count <- length(condor_clean)
     matching_cell_types_count <- length(matching_cell_types)



     # Check if there are cell types not matching
     diff_count <- total_cell_types_count - matching_cell_types_count

     if (diff_count != 0) {

       missing_cell_types <- setdiff(names(condor_clean), matching_cell_types)

       warning(paste("There are ", diff_count, " cell types that are not matching between Condor and FlowJo.\n",
                      "Missing cell types: ", paste(missing_cell_types,collapse = ", ")))

       if(diff_count == total_cell_types_count){

         stop("No matching cell types found. Ensure column names are formatted correctly and match between Condor and FlowJo.")
       }
     }


     # Subset the two data frames only for the matching cell types
     selected_condor <- condor_clean[matching_cell_types]
     selected_flowjo <- flowjo_clean[matching_cell_types]

     # Create a dataframe with only the selected columns
     cor_data <- merged_data[, c(selected_condor, selected_flowjo)]


     M <- cor(cor_data,use = "pairwise.complete.obs",method = method_corr)

     corr_M <- M[selected_condor,selected_flowjo]
     corrplot::corrplot(corr_M,method = "circle",type = "lower", tl.cex = 3.5,cl.cex = 4.5)

     plot <- recordPlot()

     return(plot)
}
