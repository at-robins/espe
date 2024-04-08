FROM python:3.9.18-bookworm
RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y autoconf automake make gcc g++ pkg-config libhdf5-dev llvm-14-dev libssl-dev libcurl4-openssl-dev less locales wget ca-certificates fonts-texgyre libharfbuzz-dev libfribidi-dev r-base gfortran libglpk-dev libgmp-dev libxml2-dev

RUN yes | pip install --upgrade pip
RUN yes | pip install numpy==1.26.1 --exists-action i
RUN yes | pip install seaborn==0.12.2 --exists-action i
RUN yes | pip install rpy2[all]==3.5.14 --exists-action i
RUN yes | pip install pathvalidate==3.2.0 --exists-action i

# Install R dependencies.
# The library calls verify that the installation worked.
WORKDIR /
ENV R_REPO="\"https://cloud.r-project.org/\""
RUN Rscript -e "install.packages(\"remotes\", repos = $R_REPO)"
RUN Rscript -e "library(remotes)"
RUN Rscript -e "remotes::install_version(\"BiocManager\", version = \"1.30.22\", repos = $R_REPO)"
RUN Rscript -e "BiocManager::install(\"DESeq2\", version = \"3.18\")"
RUN Rscript -e "library(DESeq2)"

COPY run_differential_accessibility.py /run.py
ENTRYPOINT ["python", "run.py"]