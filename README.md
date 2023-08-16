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
different data types (e.g. MCFC, CyTOF, spectral Flow). *condor* can also directly read a FlowJo workspace and apply the same biex trasfromation to the data.

**Dimensionality Reduction:** Several dimensionality reductions
algorithms are provided with the *condor* package: **PCA, UMAP, DM**.

**Batch correction:** *condor* implements *harmony* for integration of
the data.

**Clustering:** Phenograph clustering and Flow SOM clustering.

**Pseudotime:** Pseudotime can be calculated and visualized.

**Statistical Analysis:** Differential protein expression and changes in
frequencies can be tested with built in functions.

**Data Visualization:** *condor* includes several data visualization
function including: *confusion matrix, marker heatmap, feature plot,
boxplot*.

**Scalability to large scale datasets:** *condor* can work with large scale dataset, 100+M cells and can integrate new data in existing models.

**Disease classifier:** *condor* can be applied to a clinical with a classifier pre-trained.

**Misc:** *condor* includes functions to check the integrity of the flow
cytometry dataset, subset it or annotate metaclusters. Also includes several visualization tools to help in the interpretation of the data.

**On hold**
- Update to the latest version of R 4XX
- Astir for automated cell type annotation - Have some code 
- ShinyApp (student project)

## Contact or follow us

For any problem or question regarding the *condor* package or this
repository or you just want to be up to date on what is coming next
follow me:

<img src="./.logo/twitter.png" width="12%" style="float: left;"/>

[@LorenzoBonaguro](<https://twitter.com/LorenzoBonaguro>)