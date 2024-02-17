# cyCONDOR (HDC analysis ecosystem)

<img src="man/figures/condor_logo_new.png" alt="drawing" width="200" align="right"/>

High-dimensional cytometry (HDC) is a powerful tool for studying single-cells phenotypes in a complex system. Although in recent years the combination of technological developments and affordability have made HDC broadly available, these technological advances were not paired with the adequate development of analytical methods to take full advantage of the data generated. While several platforms and bioinformatics tools are currently available for the analysis of HDC data, they are either web-hosted with limited scalability or designed for expert computational biologists, making them difficult to approach by wet lab scientists. Additionally, the need for end-to-end HDC data analysis tools within a unified ecosystem poses a significant challenge, as researchers must navigate multiple platforms and software packages to complete the analysis.
We developed an easy-to-use computational framework (condor) covering not only all of the essential steps of cytometry data analysis but also including an array of downstream functions and tools to expand the biological interpretation of the data. condor's comprehensive suite of features, including guided pre-processing, clustering, dimensionality reduction, and machine learning algorithms, facilitates the seamless integration of condor into clinically relevant settings, where scalability and disease classification are paramount for the widespread adoption of HDC in clinical practice. Additionally, condor's advanced analytical features, such as pseudotime analysis and batch integration, provide researchers with the tools to extract deeper insights from their data. We used condor on a variety of data from different tissues and technologies demonstrating its versatility to assist the analysis of high dimensionality data from preprocessing to biological interpretation.

## How to use

We recommend using `cyCONDOR` from our pre-build `Docker` container [lorenzobonaguro/condor](https://hub.docker.com/r/lorenzobonaguro/condor), the latest version of the image can be pulled with:
```
docker pull lorenzobonaguro/cycondor:latest
```

To run the image you can then follow the following script

```
docker run -dp [YOURPORT]:8787 \
-e USER=[YOURUSERNAME] -e PASSWORD=[YOURPASSWORD] \
--name condor_analysis \
-v [PATHTODATA]:/home/[YOURUSERNAME]/data/ \
lorenzobonaguro/cycondor:latest
```

If you plan to run the image on `Singularity` instead of `Docker` feel free to contact us for support in the configuration.

## How to install locally

The tools was tested with `R v4.3.X`, older version should be compatible but were not tested

To install `cyCONDOR` locally first you need to install all the dependencies, from R execute the following. This script will install add dependencies from Bioconductor, CRAN and GitHub. For some package a compiler is required (e.g. Rtools on Windows or Xcode on MacOS)
```
download.file(url = "https://github.com/lorenzobonaguro/cyCONDOR/blob/master/inst/install_locally_script.R", 
              destfile = "install_locally_script.R")
              
source(install_locally_script.R)
``` 

Then install the latest release of `cyCONDOR` with:
```
devtools::install_url("https://github.com/lorenzobonaguro/cyCONDOR/releases/download/v015/cyCONDOR_0.1.5.tar.gz")
```

## Key cyCONDOR features include:

**Data loading:** Loading of `.fcs` or `.csv` files from a folder and match it with an annotation table. The data are structure in a *flow cytometry dataset* which is the structure needed for all the downstream functions. Also `FlowJo Archives` can be directly imported in cyCONDOR

**Data Transformation:** Data are transformed to be suitable for the downstream analysis, different transformation methods are provided for different data types (e.g. MCFC, CyTOF, Spectral Flow, CITE-seq).

**Dimensionality Reduction:** Several dimensionality reductions
algorithms are provided with the `cyCONDOR` package: **PCA, UMAP, tSNE, DM**.

**Batch correction:** `cyCONDOR` implements `harmony` for data integration.

**Clustering:** `Phenograph` clustering and `FlowSOM` clustering are included.

**Pseudotime:** Pseudotime and trajectory analysis with `slingshot`.

**Statistical Analysis:** Differential protein expression and changes in
frequencies can be tested with built in functions.

**Data Visualization:** `cyCONDOR` includes several data visualization
function including: *confusion matrix, marker heatmap, feature plot,
boxplot* and other visualization to also depict the data similarly to conventional cytometry data analysis software (e.g. FlowJo).

**Scalability to large scale datasets:** `cyCONDOR` can work with large scale datasets, and can integrate new data in existing models allowing population level integrated data analysis.

**Disease classifier:** `cyCONDOR` can be applied as clinical classifier with easy to train machine learning models.

**Misc:** `cyCONDOR` includes functions to check the integrity of the flow
cytometry dataset, subset it or annotate metaclusters. Also includes several visualization tools to help in the interpretation of the data. Feel free to browse all the vignette in the *Articles* section.

## How to cite cyCONDOR

> **Unveiling the Power of High-Dimensional Cytometry Data with condor**, in preparation

## Contact or follow us

For any problem or question regarding the `cyCONDOR` package or this
repository or you just want to be up to date on what is coming next
follow me:

<img src="man/figures/twitter.png" width="8%" style="float: left;"/>

[@LorenzoBonaguro](<https://twitter.com/LorenzoBonaguro>)
