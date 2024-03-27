FROM python:3.9.18-bookworm
RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y autoconf automake make gcc g++ pkg-config libhdf5-dev llvm-14-dev libssl-dev libcurl4-openssl-dev less locales wget ca-certificates fonts-texgyre libharfbuzz-dev libfribidi-dev r-base gfortran libglpk-dev libgmp-dev libxml2-dev

# Installs python dependencies.
RUN yes | pip install --upgrade pip
RUN yes | pip install uropa==4.0.3 --exists-action i

# Installs R dependencies.
# The library calls verify that the installation worked.
ENV R_REPO="\"https://cloud.r-project.org/\""
RUN Rscript -e "install.packages(\"remotes\", repos = $R_REPO)"
RUN Rscript -e "library(remotes)"
RUN Rscript -e "remotes::install_version(\"BiocManager\", version = \"1.30.22\", repos = $R_REPO)"
RUN Rscript -e "remotes::install_version(\"ggplot2\", version = \"3.5.0\", repos = $R_REPO)"
RUN Rscript -e "library(ggplot2)"
RUN Rscript -e "remotes::install_version(\"devtools\", version = \"2.4.5\", repos = $R_REPO)"
RUN Rscript -e "library(devtools)"
RUN Rscript -e "remotes::install_version(\"gplots\", version = \"3.1.3.1\", repos = $R_REPO)"
RUN Rscript -e "library(gplots)"
RUN Rscript -e "remotes::install_version(\"gridExtra\", version = \"2.3\", repos = $R_REPO)"
RUN Rscript -e "library(gridExtra)"
RUN Rscript -e "remotes::install_version(\"jsonlite\", version = \"1.8.8\", repos = $R_REPO)"
RUN Rscript -e "library(jsonlite)"
RUN Rscript -e "remotes::install_version(\"VennDiagram\", version = \"1.7.3\", repos = $R_REPO)"
RUN Rscript -e "library(VennDiagram)"
RUN Rscript -e "remotes::install_version(\"snow\", version = \"0.4-4\", repos = $R_REPO)"
RUN Rscript -e "library(snow)"
RUN Rscript -e "remotes::install_version(\"getopt\", version = \"1.20.4\", repos = $R_REPO)"
RUN Rscript -e "library(getopt)"
RUN Rscript -e "remotes::install_version(\"tidyr\", version = \"1.3.1\", repos = $R_REPO)"
RUN Rscript -e "library(tidyr)"
RUN Rscript -e "remotes::install_version(\"UpSetR\", version = \"1.4.0\", repos = $R_REPO)"
RUN Rscript -e "library(UpSetR)"
RUN Rscript -e "BiocManager::install(\"RBGL\", version = \"3.16\")"
RUN Rscript -e "library(RBGL)"
RUN Rscript -e "BiocManager::install(\"graph\", version = \"3.16\")"
RUN Rscript -e "library(graph)"
RUN Rscript -e "library(devtools) ; install_github(\"js229/Vennerable\", ref = \"46057c9\");"
RUN Rscript -e "library(Vennerable)"

COPY run_peak_annotation.py /run.py
ENTRYPOINT ["python", "run.py"]