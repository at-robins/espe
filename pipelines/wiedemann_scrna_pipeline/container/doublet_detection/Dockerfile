FROM python:3.9.18-bookworm
RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y autoconf automake make gcc g++ pkg-config libhdf5-dev llvm-14-dev libssl-dev libcurl4-openssl-dev less locales wget ca-certificates fonts-texgyre libharfbuzz-dev libfribidi-dev r-base

RUN yes | pip install --upgrade pip
RUN yes | pip install numpy==1.26.1 --exists-action i
RUN yes | pip install scanpy[leiden]==1.9.6 --exists-action i
RUN yes | pip install seaborn==0.12.2 --exists-action i
RUN yes | pip install rpy2[all]==3.5.14 --exists-action i
RUN yes | pip install anndata2ri==1.3.1 --exists-action i


# Install R dependencies.
# The library calls verify that the installation worked.
WORKDIR /
ENV R_REPO="\"https://cloud.r-project.org/\""
RUN Rscript -e "install.packages(\"remotes\", repos = $R_REPO)"
RUN Rscript -e "library(remotes)"
RUN Rscript -e "remotes::install_version(\"Matrix\", version = \"1.6-1\", repos = $R_REPO)"
RUN Rscript -e "library(Matrix)"
RUN Rscript -e "remotes::install_version(\"BiocManager\", version = \"1.30.22\", repos = $R_REPO)"
RUN Rscript -e "BiocManager::install(\"scater\", version = \"3.16\")"
RUN Rscript -e "library(scater)"
RUN Rscript -e "BiocManager::install(\"scDblFinder\", version = \"3.16\")"
RUN Rscript -e "library(scDblFinder)"
RUN Rscript -e "remotes::install_version(\"Seurat\", version = \"5.0.0\", repos = $R_REPO)"
RUN Rscript -e "install.packages(\"Seurat\", repos = $R_REPO)"
RUN Rscript -e "library(Seurat)"

COPY run_doublet_detection.py /run.py
ENTRYPOINT ["python", "run.py"]