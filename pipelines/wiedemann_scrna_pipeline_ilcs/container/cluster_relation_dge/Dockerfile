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
WORKDIR /
ENV R_REPO="\"https://cloud.r-project.org/\""
ENV BIOCONDUCTOR_VERSION="\"3.16\""
ENV R_PROJECT_PATH="\"/renv\""
RUN Rscript -e "install.packages(\"remotes\", repos = $R_REPO)" && Rscript -e "library(remotes)"
RUN Rscript -e "remotes::install_version(\"BiocManager\", version = \"1.30.23\", repos = $R_REPO)"
RUN Rscript -e "remotes::install_version(\"renv\", version = \"1.0.7\", repos = $R_REPO)" && Rscript -e "library(renv)"
RUN Rscript -e "renv::init(project = $R_PROJECT_PATH, repos = $R_REPO, bioconductor = $BIOCONDUCTOR_VERSION)"
RUN Rscript -e "renv::install(\"Matrix@1.6-5\", prompt = FALSE, repos = $R_REPO, lock = TRUE, project = $R_PROJECT_PATH)"
RUN Rscript -e "renv::install(\"SeuratObject@5.0.2\", prompt = FALSE, repos = $R_REPO, lock = TRUE, project = $R_PROJECT_PATH)"
RUN Rscript -e "renv::install(\"Seurat@5.0.2\", prompt = FALSE, repos = $R_REPO, lock = TRUE, project = $R_PROJECT_PATH)"
RUN Rscript -e "renv::install(\"bioc::scater\", prompt = FALSE, repos = $R_REPO, lock = TRUE, project = $R_PROJECT_PATH)"
RUN Rscript -e "renv::install(\"bioc::edgeR\", prompt = FALSE, repos = $R_REPO, lock = TRUE, project = $R_PROJECT_PATH)"
# Reinstalling irlba to prevent Matrix error.
# See: https://github.com/bwlewis/irlba/issues/70
RUN Rscript -e "renv::install(\"irlba@2.3.5.1\", repos = $R_REPO, lock=TRUE, project = $R_PROJECT_PATH)"


COPY run_cluster_relation_dge.py /run.py
ENTRYPOINT ["python", "run.py"]