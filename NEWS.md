# cyCONDOR 0.3.1
* Added function `condor_info` to easily explore the content of a condor object.
* Solved Issue #31, color pallette is reversed in `plot_marker_ridgeplot()` function.
* Solved Issue #30, the pseudotime function can now use a `prefix` argument.
* Improved documentation regarding the use of different co-factors when using the archsin data transformation.
* Solved Issue #33, the `plot_marker_dotplot` function have now an `alpha` argument to control the transparency of the dots.
* Solved Issue #34, the `plot_marker_ridgeplot` now have a `title` argument.
* Solved Issue #19, if an empty `fcs` or `csv` file is provided the `prep_fcd` function will stop and point to the corrupted file.
* Added function `order_param` to allow manual ordering of an annotation or clustering variable.
* Added `prefix` argument to `astir` functions.

# cyCONDOR 0.3.0
* Bug fix with data projection where specific parameters were removed in UMAP calculation
* Updated official Docker image to the latest Bioconductor version (Bioconductor 3.20, R 4.4.2)
* Completed documentation of the `Astir` workflow.
* Added parameter to select number of PCs and markers for UMAP projection in `learnUMAP` function.
* Added documentation and guide to install cyCONDOR in Silicon Mac. We now how a compiled version of the Rphenoannoy package.
* Simplified installation process, now automatically installing also the GitHub dependencies.
* Enable Geometric Sketching to subset a `condor` object with the function `geometric_sketching`
* Developed a function (`condor_session_info`) to include the `SessionInfo` to the `condor` object to increase tackability of the analysis.
* Included storage of the statistical test results in the `condor` object
* Added the option to display statistics to the `plot_frequency_boxplot` function - so far for `wilcox`, `t_test` and `diffcyt`
* Added function to correlate the cyCONDOR results and the manual gating with FlowJo or similar software (`corr_plot_comparison`)
* **Experimental**: Added a function to write FCS files from a `condor` object (`write_fcs`). Currently the function if writing in the .fcs file the cyCONDOR transformed values, we will change to the original values in the next iteration after more extensive testing.

# cyCONDOR 0.2.1
* Included help function to assign metaclusters (Thanks to Lucas Secchim Ribeiro)
* Removal of redundant `color_palette` parameter and other redundant lines in the source code
* Added `turkey hsd` as post-hoc test to `frequency_anova_test` (see issue #10)
* Modified `emmeans test` in `frequency_anova_test` (see issue #10)
* Improvement of error messages

# cyCONDOR 0.2.0
* Bump of major version for official release and publication
* Changed order of events in the `prep_fcd` function. This makes the process faster (especially in smaller machines) and reduces the usage of memory.
* Added function `subset_fcd_byparam`. This function allows to randomly subset a `condor` object proportionally across a selected parameter. For example is not possible to randomly subset `n` cells from each samples.
* Included more professionalization in the `read_fcd` function. This allows the user to customize few aspects of the `condor` object
* Included `cyCONDOR` version in the `condor` object under the `extras` slot. This enables the user to trace the version used for the analysis.

# cyCONDOR 0.1.6

* Reorganization of existing visualization functions including harmonization of function names and function arguments, utilization of `condor` object as main input object and addition of more extensive documentation and error messages.
* Added visualization functions `plot_counts_barplot()`, `plot_marker_ridgeplot()` and `plot_marker_boxplot()`
* Added `getTable()` function to generate tables of cell population counts and frequencies, as well as mean or median marker expression for all cell population - sample - marker combinations.
* `boxplot_and_stats()` function was replaces by `plot_frequency_boxplot()` function for visualization and several functions to conduct statistical tests on population frequencies.
* Added wrapper functions around basic statistical tests to compare cell population frequencies between groups of samples (`frequency_t_test()`, `frequency_wilcox_test()`, `frequency_anova_test()`, `frequency_kruskal_test()`, `frequency_friedman_test()`) 
* Added `prepInputDiffcyt()` function to transform the `condor` object into an SummarizedExperiment object compatible with the `diffcyt` package for differential testing.
* Renaming of arguments in `runPseudotime()` function to harmonize within the package
* Renaming of arguments in `metaclustering()` function to harmonize within the package
* Updated documentation
* Renaming of arguments in multiple function to harmonize within the package
* Setting a default seed in multiple functions
* Added functions to use the `CytoNorm` algorithm for batch normalization
* Bug fixes in `prep_flw()` when merging annotation and removing parameters, saving of import parameters in extras slot of fcd
* Added new parameters in `runFlowSOM()` to determine size of FlowSOM grid
* implemented marker selection (`runPCA()`, `runUMAP()`, `runDM()`, `runtSNE()`, `runPhenograph()`, `runFlowSOM()`), saving of marker selection in extra slot of fcd
* Added functions to extract all markers present in fcd (`measured_markers()`) or selected markers (`used_markers()`)
* Added function to visualize PC loadings
* Simplified data loading and transformation function including useful error messages
* Included `arcsinh` transformation with cofactor 5 for cyTOF data

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
