FROM python:3.11.12-bookworm
RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y autoconf automake make gcc g++ pkg-config libhdf5-dev llvm-14-dev libssl-dev libcurl4-openssl-dev less locales wget ca-certificates fonts-texgyre libharfbuzz-dev libfribidi-dev r-base gfortran libglpk-dev libgmp-dev libxml2-dev

RUN yes | pip install --upgrade pip
RUN yes | pip install numpy==2.2.5 --exists-action i
RUN yes | pip install pandas==2.2.3 --exists-action i
RUN yes | pip install seaborn==0.13.2 --exists-action i
RUN yes | pip install pathvalidate==3.2.3 --exists-action i
RUN yes | pip install scipy==1.15.2 --exists-action i
RUN yes | pip install anndata==0.11.4 --exists-action i
RUN yes | pip install liana==1.5.1 --exists-action i
RUN yes | pip install git+https://github.com/saezlab/decoupler-py.git@2d97793bf2b35b9926576be28028a591b87ff03a

# Downloads the required pathway network files.
COPY download_pathways.py /download_pathways.py
RUN python /download_pathways.py

COPY run_pathway_analysis.py /run.py
ENTRYPOINT ["python", "/run.py"]