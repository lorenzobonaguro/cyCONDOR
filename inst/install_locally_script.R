# Install cyCONDOR locally

###################################################################################################################
# NOTE                                                                                                            #
#                                                                                                                 #
# To install cyCONDOR on Windows or Mac you need to have a compiler installed to build the packages from source.  #
# On Windows install Rtools, the version should match the version of the major R release used                     #
# On Mac install xcode                                                                                            #
# This code was tested on R 4.4 and Bioconductor 3.20                                                             #
###################################################################################################################

# First we make sure Bioconductor is installed and updated
BiocManager::install(update = T, ask = F, version = "3.20")

# Install devtools
BiocManager::install("devtools",version = "3.20")

# Next we install two dependencies which are only available on GitHub
devtools::install_github(repo = c("JinmiaoChenLab/Rphenograph", "stuchly/Rphenoannoy", "saeyslab/CytoNorm@362ac08"),
                         repos = BiocManager::repositories())

# Finally we install cyCONDOR, here we manually provide the link to the Bioconductor repositories.
devtools::install_url("https://github.com/lorenzobonaguro/cyCONDOR/releases/download/v021/cyCONDOR_0.2.2.tar.gz",
                      repos = BiocManager::repositories())
