FROM python:3.9.18-bookworm
RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y autoconf automake make gcc g++ pkg-config libhdf5-dev llvm-14-dev libssl-dev libcurl4-openssl-dev less locales wget ca-certificates fonts-texgyre libharfbuzz-dev libfribidi-dev r-base gfortran libglpk-dev libgmp-dev libxml2-dev

RUN yes | pip install --upgrade pip
RUN yes | pip install numpy==1.26.1 --exists-action i
RUN yes | pip install scanpy[leiden]==1.9.6 --exists-action i
RUN yes | pip install seaborn==0.12.2 --exists-action i
RUN yes | pip install pathvalidate==3.2.0 --exists-action i
RUN yes | pip install rpy2[all]==3.5.14 --exists-action i
RUN yes | pip install anndata2ri==1.3.1 --exists-action i

# Install R dependencies.
# The library calls verify that the installation worked.
ENV R_REPO="\"https://cloud.r-project.org/\""
ENV BIOCONDUCTOR_VERSION="\"3.16\""
ENV R_PROJECT_PATH="/renv/"
RUN Rscript -e "install.packages(\"remotes\", repos = $R_REPO)" && Rscript -e "library(remotes)"
RUN Rscript -e "remotes::install_version(\"BiocManager\", version = \"1.30.23\", repos = $R_REPO)"
RUN Rscript -e "remotes::install_version(\"renv\", version = \"1.0.7\", repos = $R_REPO)" && Rscript -e "library(renv)"
COPY r_sc_find_pca.lock $R_PROJECT_PATH
RUN Rscript -e "renv::restore(\
    lockfile = \"${R_PROJECT_PATH}r_sc_find_pca.lock\", \
    library = \"/usr/local/lib/R/site-library\", \
    repos = $R_REPO, \
    prompt = FALSE \
)"
RUN Rscript -e "library(findPC)"

# This section can be used to install R dependencies and output the lock file.
# RUN Rscript -e "install.packages(\"remotes\", repos = $R_REPO)" && Rscript -e "library(remotes)"
# RUN Rscript -e "remotes::install_version(\"BiocManager\", version = \"1.30.23\", repos = $R_REPO)"
# RUN Rscript -e "remotes::install_version(\"renv\", version = \"1.0.7\", repos = $R_REPO)" && Rscript -e "library(renv)"
# RUN Rscript -e "renv::init(project = \"$R_PROJECT_PATH\", repos = $R_REPO, bioconductor = $BIOCONDUCTOR_VERSION)"
# RUN Rscript -e "renv::install(\"haotian-zhuang/findPC@0e85c94\", prompt = FALSE, repos = $R_REPO, lock = TRUE, project = \"$R_PROJECT_PATH\")"
# RUN cat /renv/renv.lock

COPY run_clustering_cell_type.py /run.py
ENTRYPOINT ["python", "run.py"]