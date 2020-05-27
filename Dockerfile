# default R-version install via rocker, uses debian as a base
FROM rocker/r-ver:3.5.2

# install some likely used base unix dependencies
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
	libxml2-dev \
	libxt-dev \
	libjpeg-dev \
	libglu1-mesa-dev \
	libcairo2-dev \
	libsqlite3-dev \
	libmariadbd-dev \
	libmariadb-client-lgpl-dev \
	libpq-dev \
	libmagick++-dev \
	libssh2-1-dev \
	libssl-dev \
	libcurl4-openssl-dev \
	libnss3 \
	libclang-dev \
	unixodbc-dev \
	cargo \
	wget \
	bzip2 \
	nano
	
# create some basic directories for binding ICB storages
RUN mkdir -p /storage/groups/ /storage/scratch/ /home/icb/johann.hawe/

# install some custom software
RUN mkdir /software/ && cd /software/

# install some basic packages
RUN install2.r --error \
	--deps TRUE \
	plyr \
	dplyr \
	knitr \
	tidyverse \
	DT \
	data.table \
	Rcpp 

# also add some packages which just could be really useful
RUN R -e "source('https://bioconductor.org/biocLite.R')" \
	-e "biocLite(c('feather', 'ggpubr', 'GenomicRanges', 'batchtools', 'mclust', 'mixR', 'mixtools'))"
	
RUN install2.r --error \
	--deps TRUE \
	devtools

# add custom packages
RUN R -e "source('https://bioconductor.org/biocLite.R')" \
	-e "biocLite(c('BDgraph', 'GENIE3', 'GeneNet', 'Homo.sapiens', 'httr', 'Matrix', 'org.Hs.eg.db'))"

RUN R -e "source('https://bioconductor.org/biocLite.R')" \
	-e "biocLite(c('pheatmap', 'preprocessCore', 'plsgenomics', 'qvalue', 'fdrtool', 'reshape2'))"

RUN R -e "source('https://bioconductor.org/biocLite.R')" \
	-e "biocLite(c('RBGL', 'rtracklayer', 'scales', 'sva', 'illuminaHumanv3.db', 'FDb.InfiniumMethylation.hg19'))"

RUN R -e "source('https://bioconductor.org/biocLite.R')" \
	-e "biocLite(c('AnnotationDbi', 'annotables', 'BSgenome.Hsapiens.UCSC.hg19'))"

RUN R -e "source('https://bioconductor.org/biocLite.R')" \
	-e "biocLite(c('Rsamtools', 'Rgraphviz', 'cvTools', 'ROCR', 'doRNG'))"

RUN R -e "source('https://bioconductor.org/biocLite.R')" \
	-e "biocLite(c('GSEABase', 'GOstats'))"

RUN R -e "source('https://bioconductor.org/biocLite.R')" \
	-e "biocLite(c('glasso'))"
	
RUN mkdir /r-packages/
COPY iRafNet_1.1-2.tar.gz /r-packages/iRafNet_1.1-2.tar.gz
RUN R CMD INSTALL /r-packages/iRafNet_1.1-2.tar.gz

RUN R -e "source('https://bioconductor.org/biocLite.R')" \
	-e "biocLite(c('huge'))"

RUN R -e "library(devtools)" \
	-e "install_github('zdk123/SpiecEasi')"

RUN R -e "library(devtools)" \
    -e "devtools::install_github('andreyshabalin/MatrixEQTL')"

RUN R -e "source('https://bioconductor.org/biocLite.R')" \
	-e "biocLite(c('meta'))"
	
RUN R -e "source('https://bioconductor.org/biocLite.R')" \
	-e "biocLite(c('mvtnorm'))"

# system clean up
RUN rm -rf /var/lib/apt/lists/* \
	&& apt-get clean \
	&& apt-get purge
	
# R cleanup
RUN rm -rf /tmp/downloaded_packages/ /tmp/*.rds
