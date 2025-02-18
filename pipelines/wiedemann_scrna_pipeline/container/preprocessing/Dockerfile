FROM python:3.9.18-bookworm
RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y autoconf automake make gcc g++ pkg-config libhdf5-dev llvm-14-dev libssl-dev libcurl4-openssl-dev less locales wget ca-certificates fonts-texgyre libharfbuzz-dev libfribidi-dev r-base gfortran libglpk-dev libgmp-dev libxml2-dev

RUN yes | pip install --upgrade pip
RUN yes | pip install numpy==1.26.4 --exists-action i
RUN yes | pip install pandas==2.2.2 --exists-action i
RUN yes | pip install seaborn==0.13.2 --exists-action i
RUN yes | pip install pathvalidate==3.2.0 --exists-action i
RUN yes | pip install scanpy[leiden]==1.9.6 --exists-action i
RUN yes | pip install Cython==3.0.12 --exists-action i
RUN yes | pip install pegasuspy==1.10.1 --exists-action i

RUN git clone https://github.com/ayshwaryas/ddqc.git && git -C /ddqc checkout b09a71e && pip install /ddqc

COPY run_preprocessing.py /run.py
ENTRYPOINT ["python", "run.py"]