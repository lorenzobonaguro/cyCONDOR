# This docker image is meant for the analysis of high dimensionality flow data with CONDOR

# I start from the Bioconductor 3.17 image with RStudio Server
FROM bioconductor/bioconductor_docker:RELEASE_3_20

# Install miniconda
ENV CONDA_DIR /opt/conda
RUN wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh" -O ~/miniforge.sh && \
    /bin/bash ~/miniforge.sh -b -p /opt/conda

RUN rm /root/miniforge.sh

# Put conda in path so we can use conda activate
ENV PATH=$CONDA_DIR/bin:$PATH

# Add the conda repositories
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda
RUN conda config --set channel_priority strict

## Install astir in this conda
RUN cd /opt && git clone https://github.com/camlab-bioml/astir.git && cd astir && conda run -n astir pip install .

RUN pip install geosketch

# Update Bioconductor and install devtools
RUN Rscript -e 'BiocManager::install(update = T, ask = F, version = "3.20")'

RUN Rscript -e 'BiocManager::install("devtools",version = "3.20")'

# Install the GitHib dependencies
RUN Rscript -e 'devtools::install_github(repo = c("JinmiaoChenLab/Rphenograph", "stuchly/Rphenoannoy", "saeyslab/CytoNorm@362ac08"), repos = BiocManager::repositories())'

# install CONDOR
RUN Rscript -e 'devtools::install_url("https://github.com/lorenzobonaguro/cyCONDOR/releases/download/v021/cyCONDOR_0.2.1.tar.gz", repos = BiocManager::repositories())'

# example to run this docker image
# docker run -dp 8787:8787 -e PASSWORD=mariorossi --name fc_analysis -v /mnt/e/docker_test:/home/lorenzo/data/ lorenzobonaguro/condor:v022
