FROM python:3.11.12-bookworm
RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y autoconf automake make gcc g++ pkg-config libhdf5-dev llvm-14-dev libssl-dev libcurl4-openssl-dev less locales wget ca-certificates fonts-texgyre libharfbuzz-dev libfribidi-dev r-base gfortran libglpk-dev libgmp-dev libxml2-dev

RUN yes | pip install --upgrade pip
RUN yes | pip install numpy==2.2.5 --exists-action i
RUN yes | pip install beautifulsoup4==4.13.4 --exists-action i
RUN yes | pip install khmer==2.1.1 --exists-action i

COPY run_mappable_size.py /run.py
ENTRYPOINT ["python", "/run.py"]