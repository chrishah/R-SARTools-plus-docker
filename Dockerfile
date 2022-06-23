FROM rocker/r-ver:4.0.3

RUN apt-get update && \
	apt-get install -y \
		zlib1g-dev \
		libxml2 \
		libbz2-dev \
		liblzma-dev \
		libpcre3-dev \
		libicu-dev \
		libcairo2 \
		libcairo2-dev \
		libxt-dev 


#install r packages
RUN install2.r --error \
	devtools \
	XML \
	dplyr \
	BiocManager \
	Cairo

#install R packages via devtools
RUN R -e "devtools::install_github('PF2-pasteur-fr/SARTools@2b95eaa473c20d85f0fa92305acdf43f68f9baed', build_opts='--no-resave-data')" && \
	R -e "BiocManager::install('ComplexHeatmap')"

RUN apt-get install -y pandoc texlive
RUN apt-get install -y pandoc-citeproc
RUN install2.r --error gplots
