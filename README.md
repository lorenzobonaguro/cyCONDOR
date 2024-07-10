# cyCONDOR (HDC analysis ecosystem)

<img src="man/figures/condor_logo_new.png" alt="drawing" width="200" align="right"/>

High-dimensional cytometry (HDC) is a powerful tool for studying single-cells phenotypes in a complex system. Although in recent years the combination of technological developments and affordability have made HDC broadly available, these technological advances were not paired with the adequate development of analytical methods to take full advantage of the data generated. While several platforms and bioinformatics tools are currently available for the analysis of HDC data, they are either web-hosted with limited scalability or designed for expert computational biologists, making them difficult to approach by wet lab scientists. Additionally, the need for end-to-end HDC data analysis tools within a unified ecosystem poses a significant challenge, as researchers must navigate multiple platforms and software packages to complete the analysis.
We developed an easy-to-use computational framework (condor) covering not only all of the essential steps of cytometry data analysis but also including an array of downstream functions and tools to expand the biological interpretation of the data. condor's comprehensive suite of features, including guided pre-processing, clustering, dimensionality reduction, and machine learning algorithms, facilitates the seamless integration of condor into clinically relevant settings, where scalability and disease classification are paramount for the widespread adoption of HDC in clinical practice. Additionally, condor's advanced analytical features, such as pseudotime analysis and batch integration, provide researchers with the tools to extract deeper insights from their data. We used condor on a variety of data from different tissues and technologies demonstrating its versatility to assist the analysis of high dimensionality data from preprocessing to biological interpretation.

## How to use

**You can find the detailed documentation of cyCONDOR [here](https://lorenzobonaguro.github.io/cyCONDOR/)**

We recommend using `cyCONDOR` from our pre-build `Docker` container [lorenzobonaguro/condor](https://hub.docker.com/r/lorenzobonaguro/condor), the latest version of the image can be pulled with:
```
docker pull lorenzobonaguro/cycondor:v016
```

To run the image you can then follow the following script

```
docker run -dp [YOURPORT]:8787 \
-e USER=[YOURUSERNAME] -e PASSWORD=[YOURPASSWORD] \
--name condor_analysis \
-v [PATHTODATA]:/home/[YOURUSERNAME]/data/ \
lorenzobonaguro/cycondor:v016
```
You can then access RStudio from your web browser at the address

```
http://localhost:[YOURPORT]/
```

If you are starting the Docker container from a remote server exchange the `localhost` with the IP address or domain name of the server as exemplified below:

```
http://[SERVERNAME]:[YOURPORT]/
```

A detailed guide on how to get started with Docker and how to run cyCONDOR as `Singulariy` container are provided in the vignette: `vignette("How_to_run_cyCONDOR_as_container")`.

## How to install locally

The tools was tested with `R v4.3.X`, older version should be compatible but were not tested

To install `cyCONDOR` you can follow few steps describe here below. 

**IMPORTANT:** For some package a compiler is required (e.g. Rtools on Windows or Xcode on MacOS)

First install `Bioconductor`, if you are sure `Bioconductor` is already installed in your system you can skip this step.
```
BiocManager::install(update = T, ask = F, version = "3.17")
``` 

Next we install two dependencies which are only available on GitHub
```
devtools::install_github(repo = c("JinmiaoChenLab/Rphenograph", "stuchly/Rphenoannoy", "saeyslab/CytoNorm@362ac08"))
```

Finally we install cyCONDOR, here we manually provide the link to the Bioconductor repositories.
```
devtools::install_url("https://github.com/lorenzobonaguro/cyCONDOR/releases/download/v016/cyCONDOR_0.1.6.tar.gz",
                      repos = BiocManager::repositories())
```

Alternatively those steps could be automated with the following code
```
download.file(url = "https://raw.githubusercontent.com/lorenzobonaguro/cyCONDOR/master/inst/install_locally_script.R", destfile = "install_locally_script.R")
              
source(install_locally_script.R)
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

> [**Kroeger, Mueller, Leidner et al. Unveiling the Power of High-Dimensional Cytometry Data with cyCONDOR**](https://www.biorxiv.org/content/10.1101/2024.02.29.582727v1)

## Contact or follow us

For any problem or question regarding the `cyCONDOR` package or this
repository or you just want to be up to date on what is coming next
follow me:

<img src="man/figures/twitter.png" width="8%" style="float: left;"/>

[@LorenzoBonaguro](<https://twitter.com/LorenzoBonaguro>)
