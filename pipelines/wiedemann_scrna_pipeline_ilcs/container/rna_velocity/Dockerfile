FROM python:3.8.19-bookworm
RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y autoconf automake make gcc g++ pkg-config libhdf5-dev llvm-14-dev libssl-dev libcurl4-openssl-dev less locales wget ca-certificates fonts-texgyre libharfbuzz-dev libfribidi-dev r-base gfortran libglpk-dev libgmp-dev libxml2-dev

RUN yes | pip install --upgrade pip
RUN yes | pip install pandas==1.2.3 --exists-action i
RUN yes | pip install anndata==0.8.0 --exists-action i
RUN yes | pip install scanpy==1.8.2 --exists-action i
RUN yes | pip install numpy==1.21.6 --exists-action i
RUN yes | pip install scipy==1.10.1 --exists-action i
RUN yes | pip install numba==0.57.0 --exists-action i
RUN yes | pip install matplotlib==3.3.4 --exists-action i
RUN yes | pip install scvelo==0.2.4 --exists-action i
RUN yes | pip install typing-extensions==4.11.0 --exists-action i
RUN yes | pip install pathvalidate==3.2.0 --exists-action i

RUN git clone https://github.com/at-robins/TFvelo.git

COPY run_rna_velocity.py /run.py
ENTRYPOINT ["python", "run.py"]