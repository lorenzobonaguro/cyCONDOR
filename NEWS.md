# cyCONDOR 0.1.5

* Fix bug in the definition of tab separator when loading csv files
* Added clr transformation for CITE-seq data together with minor improvements to the transformation function
* Add clustering option in confusion matrix
* Added visualization of 2D plots of PCA
* Option to export plots as raster
* Added function to plot conventional flow 2d plots
* Added function to plot a density plot for marker expression
* Implementation of Astir (with Python)
* Restructured documentation and vignette
* Added GitPages website with documentation and tutorials
* Name change to `cyCONDOR`
* Edited the pseudotime function to run in a loop for the starting clusters

# cyCONDOR 0.1.4

## Reference name: condor 0.1.4

* Add ML classifier with CytoML
* Included Hmisc as requirment for Violin plto marker function
* Fix bug with UMAP plotting function when faceting (default for facet_by_fariable set to FALSE not NULL)
* Fixed package loading message
* Tested diffusion map and imporoved function
* Tested pseutodime and improved function
* Splittied functions in multiple files to make them easier to find
* Improved package documentation
* Improved package vignette

# cyCONDOR 0.1.3

## Reference name: condor 0.1.3

* FlowSOM function can retain the model to plot the SOM tree afterwards
* Added Function to read flowjo workspaces to a condor campatible format
* Fixed filter function when the 'extra' slot is occupied

# cyCONDOR 0.1.2

## Reference name: condor 0.1.2

* Several bug fixes
* Added function to calculate tSNE dimensionality reduction
* Possibility to limit the number of PC used for clustering and non-linear dimensionality reduction
* Added function to calculate Pseudotime (slignshot)
* Included an easy-to-export differential frequency table
* Added a function for random subsetting of the dataset
* Added visualization of PC loadings
* Added Pseudobulk PCA Analysis
* Added function to easily export cellular frequency
* Added function to change the parameter names (of the fcs files)
* Added function to visualize DRs as density plot
* Added Violin plot visualization of marker expression
* Included option to not cluster rows and columns in heatmaps
* Included option to show cluster numbers in the dotplot
* Add multicore support to tSNE
* Added workflow for UMAP projection and label transfer

# cyCONDOR 0.1.1 

## Reference name: condor 0.1.1

* Add function to merge condor objects
* Updated LoadFCS function to be fully compatible with .csv files
* Updated UMAP function to run on multiple cores
* Updated RPhenograph function to run on multiple cores
* Added function to run FlowSOM clustering
* Fix issues with UMAP parametes selections
* Added option to specifiy the delimiter for the csv files

# cyCONDOR 0.1.0

## Reference name: condor 0.1.0

* Initial release
