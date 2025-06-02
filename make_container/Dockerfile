# Use the official RStudio base image
FROM rocker/rstudio:4.2.2

# Set system locale
ENV PATH=/usr/lib/rstudio-server/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
ENV TERM=vt100
ENV LANG=en_US.UTF-8
ENV LC_ALL=en_US

# Install system dependencies
RUN apt-get -y update && apt-get install -y -qq --no-install-recommends \
        build-essential htop vim wget screen \
        libudunits2-dev libgdal-dev libgeos-dev libcairo2-dev libbz2-dev libproj-dev libxt-dev libpng-dev \
        libhdf5-dev libmagick++-dev libboost-dev libboost-all-dev libfftw3-3 libglpk40 libgsl27 libtcl8.6 && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -q -e 'install.packages(c("devtools","R.utils","Rfast","RhpcBLASctl","BiocManager"), quiet=TRUE)' && \
    R -q -e 'BiocManager::install(c("mutoss"), quiet=TRUE)' && \
    R -q -e 'install.packages(c("data.table","Matrix","matrixStats","igraph","metap","pdist"), quiet=TRUE)' && \
    R -q -e 'install.packages(c("stringr","flexmix","mixtools","qlcMatrix","hash","tm"), quiet=TRUE)' && \
    R -q -e 'install.packages(c("stringdist","aricode","bezier","dbscan","philentropy"), quiet=TRUE)' && \
    R -q -e 'install.packages(c("tidyverse","dplyr","magrittr","ggplot2","ggraph","ggrepel"), quiet=TRUE)' && \
    R -q -e 'install.packages(c("ggridges","ggbeeswarm","ggalluvial","gganimate","gifski","av"), quiet=TRUE)' && \
    R -q -e 'install.packages(c("ggdendro","patchwork","tidygraph","gridExtra","plot3D","plotly"), quiet=TRUE)' && \
    R -q -e 'install.packages(c("pheatmap","viridis","eulerr","VennDiagram","wordcloud"), quiet=TRUE)' && \
    R -q -e 'BiocManager::install(c("ComplexHeatmap"), quiet=TRUE)' && \
    R -q -e 'install.packages(c("Seurat","clustree","harmony","phateR","hdf5r","SoupX"), quiet=TRUE)' && \
    R -q -e 'devtools::install_github("mojaveazure/seurat-disk", quiet=TRUE)' && \
    R -q -e 'devtools::install_github("satijalab/seurat-wrappers", quiet=TRUE)' && \
    R -q -e 'BiocManager::install(c("SingleCellExperiment","DropletUtils","batchelor","scran"), quiet=TRUE)' && \
    R -q -e 'BiocManager::install(c("slingshot","tradeSeq","miloR","GENIE3","HCAData","sva"), quiet=TRUE)' && \
    R -q -e 'devtools::install_github("limfuxing/ruvIIInb", build_vignettes=FALSE, quiet=TRUE)' && \
    R -q -e 'BiocManager::install(c("pcaMethods","escape","UCell"), quiet=TRUE)' && \
    R -q -e 'devtools::install_github("velocyto-team/velocyto.R", quiet=TRUE)' && \
    R -q -e 'BiocManager::install("tanaylab/metacell", quiet=TRUE)' && \
    R -q -e 'BiocManager::install(c("BiocGenerics","DelayedArray","DelayedMatrixStats","limma"), quiet=TRUE)' && \
    R -q -e 'BiocManager::install(c("lme4","S4Vectors","HDF5Array","terra","ggrastr","celda"), quiet=TRUE)' && \
    R -q -e 'devtools::install_github("cole-trapnell-lab/monocle3", quiet=TRUE)' && \
    R -q -e 'devtools::install_github("saeyslab/nichenetr", quiet=TRUE)' && \
    R -q -e 'devtools::install_github("sqjin/CellChat", quiet=TRUE)' && \
    wget --no-check-certificate https://cytotrace.stanford.edu/CytoTRACE_0.3.3.tar.gz \
        -O /tmp/CytoTRACE_0.3.3.tar.gz && \
    R -q -e 'devtools::install_local("/tmp/CytoTRACE_0.3.3.tar.gz", quiet=TRUE)' && \
    R -q -e 'BiocManager::install(c("preprocessCore","DESeq2","clusterProfiler","msigdbr"), quiet=TRUE)' && \
    R -q -e 'BiocManager::install(c("fgsea","biomaRt","RUVSeq","dorothea","GenomicFeatures"), quiet=TRUE)' && \
    R -q -e 'BiocManager::install(c("Gviz","rtracklayer","BRETIGEA","sangerseqR"), quiet=TRUE)' && \
    R -q -e 'devtools::install_github("morris-lab/CellTagR", quiet=TRUE)' && \
    R -q -e 'BiocManager::install(c("cTRAP","cmapR","liftOver","affy","annotatr","genomation"), quiet=TRUE)' && \
    R -q -e 'install.packages(c("h2o","randomForest","ranger","gbm","xgboost","lightgbm","iml"), quiet=TRUE)' && \
    R -q -e 'install.packages(c("dendextend","parallelDist","metR","coop","fastcluster"), quiet=TRUE)' && \
    R -q -e 'install.packages(c("shinyhelper","DT","readr"), quiet=TRUE)' && \
    R -q -e 'devtools::install_github("SGDDNB/ShinyCell", quiet=TRUE)'

# Expose RStudio server port
EXPOSE 8787

# Start RStudio server
CMD ["rstudio-server", "start"]