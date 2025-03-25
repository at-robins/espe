FROM python:3.9.18-bookworm
RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y autoconf automake make gcc g++ pkg-config libhdf5-dev llvm-14-dev libssl-dev libcurl4-openssl-dev less locales wget ca-certificates fonts-texgyre libharfbuzz-dev libfribidi-dev r-base gfortran libglpk-dev libgmp-dev libxml2-dev

# Installs HOMER.
RUN mkdir /HOMER
WORKDIR /HOMER
RUN wget --no-verbose http://homer.ucsd.edu/homer/configureHomer.pl
RUN perl /HOMER/configureHomer.pl -install -version v5.1
RUN perl /HOMER/configureHomer.pl -install mm10 -version v6.4
RUN perl /HOMER/configureHomer.pl -install hg38 -version v6.4
ENV PATH="/HOMER/bin:${PATH}"

WORKDIR /
COPY run_motif_discovery.py /run.py
ENTRYPOINT ["python", "run.py"]