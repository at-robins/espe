FROM python:3.12.12-trixie
RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y autoconf automake make gcc g++ pkg-config libhdf5-dev llvm-19-dev libssl-dev libcurl4-openssl-dev less locales wget ca-certificates fonts-texgyre libharfbuzz-dev libfribidi-dev r-base gfortran libglpk-dev libgmp-dev libxml2-dev

RUN yes | pip install --upgrade pip
RUN yes | pip install numpy==2.4.2 --exists-action i
RUN yes | pip install scanpy[leiden]==1.12 --exists-action i
RUN yes | pip install seaborn==0.13.2 --exists-action i
RUN yes | pip install pathvalidate==3.3.1 --exists-action i
RUN yes | pip install decoupler==2.1.4 --exists-action i

COPY run_celltype_ulm.py /run.py
ENTRYPOINT ["python", "run.py"]