# WELCOME TO CONDOR (R package)

<img src="./.logo/condor_logo_new.png" alt="drawing" width="200" align="right"/>

## Overview

Flow cytometry analysis workflow. The aim of this project is to provide
an intuitive workflow for the analysis of high-dimensionality cytometry
data in R. This project is far from complete, any suggestion is well
accepted!

The tools is running on R v3.6.2, newer version are not yet fully
supported due to some problems with harmony in R v4.

An easy way to use *condor* is with the pre-built docker image
[lorenzobonaguro/condor](https://hub.docker.com/r/lorenzobonaguro/condor)

## Features currently included:

**Data loading:** Loading of .fcs files from a folder and match it with
an annotation table. The data are structure in a *flow cytometry
dataset* which is the structure needed for all the downstream functions.

**Data Transformation:** Data are transformed to be suitable for the
downstream analysis, different transformation methods are provided for
different data types.

**Dimensionality Reduction:** Several dimensionality reductions
algorithms are provided with the *condor* package: **PCA, UMAP, DM**.

**Batch correction:** *condor* implements *harmony* for integration of
the data.

**Clustering:** Phenograph clustering.

**Statistical Analysis:** Differential expression and changes in
frequencies can be tested with built in functions.

**Data Visualization:** *condor* includes several data visualization
function including: *confusion matrix, marker heatmap, feature plot,
boxplot*.

**Misc:** *condor* includes functions to check the integrity of the flow
cytometry dataset, subset it or annotate metaclusters.


## Coming soon:

**Data analysis** 
- Pseudo-time Analysis - Have some code
- Differential frequencies (boxplot of frequencies for more than two groups)
- Differential marker expression
- ~~Correct the FlowSOM~~
- ~~Function for random subsetting~~
- ~~PC loadings and variance explained by PC - Charlotte has some code~~
- Possibility to select the number of PC for clustering or DR or the features 
- ~~Pseudobulk PCA analysis~~
- Function to change param names
- Function specific to export cell frequencies

**Data Visualization** 
- ~~Density Map~~
- ~~Violin plot marker visualization~~
- ~~Cluster marker heatmap: option without column clustering~~
- Option to show cluster numbers on UMAP (Charlotte)
- ~~Violin plots of expression per cluster (sometimes I think this visualization is more helpful than the heatmap) - Charlotte has some code~~
- ~~Split UMAP by one variable and color by another~~ - Done but need further testing

**On hold**
- Import of biexponential values to recreate FlowJo visualization
- I personally find the way how names are added to the cluster numbers in the metaclustering function not so nice. I would prefer a way, where you have the number next to the label in your code.
- Update to the latest version of R 4XX
- Astir for automated cell type annotation - Have some code 

## Contact or follow us

For any problem of question regarding the *condor* package or this
repository or you just want to be up to date on what is coming next
follow me:

<img src="./.logo/twitter.png" width="12%" style="float: left;"/>

[@LorenzoBonaguro](<https://twitter.com/LorenzoBonaguro>)