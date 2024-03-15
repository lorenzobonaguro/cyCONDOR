# Install condor locally

# *NOTE*
# To install condor on Windows or Mac you need to have a compiler installed to build the packages from source.
# On Windows install Rtools, the version should match the version of the majour R release used
# On Mac install xcode
# This code was tested on R 4.3.0 and Bioconductor 3.17, nevetheless there are no expected problems with R 4.3.1 or 4.3.2 and Bioconductor 3.18

BiocManager::install(update = T, ask = F, version = "3.17")

BiocManager::install(c("flowCore", "devtools", "Rtsne", "umap", "ggplot2", "ggpubr",
                       "vegan", "diffusionMap", "FlowSOM", "gplots", "grid", "gridExtra",
                       "reshape2", "RColorBrewer", "scatterplot3d", "rgl", "Biobase", "Rmisc",
                       "tidyverse", "pheatmap", "viridis", "rainbow", "destiny",
                       "SingleCellExperiment", "slingshot", "ggrastr",
                       "uwot", "FlowSOM", "caret", "factoextra", "harmony", "flowWorkspace",
                       "CytoML", "Biobase", "Hmisc", "CytoDx", "DelayedMatrixStats"), version = "3.17")

install.packages("https://cran.r-project.org/src/contrib/randomForest_4.7-1.1.tar.gz", repos=NULL, type="source")

devtools::install_github("JinmiaoChenLab/Rphenograph")

devtools::install_github("stuchly/Rphenoannoy")

devtools::install_url("https://github.com/lorenzobonaguro/cyCONDOR/releases/download/v015/cyCONDOR_0.1.5.tar.gz")