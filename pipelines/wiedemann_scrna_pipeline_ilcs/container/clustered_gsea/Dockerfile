FROM python:3.9.18-bookworm
RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y autoconf automake make gcc g++ pkg-config libhdf5-dev llvm-14-dev libssl-dev libcurl4-openssl-dev less locales wget ca-certificates fonts-texgyre libharfbuzz-dev libfribidi-dev r-base gfortran libglpk-dev libgmp-dev libxml2-dev

RUN yes | pip install --upgrade pip
RUN yes | pip install numpy==1.26.1 --exists-action i
RUN yes | pip install scanpy[leiden]==1.9.6 --exists-action i
RUN yes | pip install seaborn==0.12.2 --exists-action i
RUN yes | pip install pathvalidate==3.2.0 --exists-action i
RUN yes | pip install decoupler==1.6.0 --exists-action i

# Downloads the required GMT files.
RUN mkdir /gmts
RUN wget -O /gmts/reactome_human.gmt https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/c2.cp.reactome.v2023.2.Hs.symbols.gmt
RUN wget -O /gmts/reactome_mouse.gmt https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Mm/m2.cp.reactome.v2023.2.Mm.symbols.gmt
RUN wget -O /gmts/go_human.gmt https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/c5.all.v2023.2.Hs.symbols.gmt
RUN wget -O /gmts/go_mouse.gmt https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Mm/m5.all.v2023.2.Mm.symbols.gmt

COPY run_clustered_gsea.py /run.py
ENTRYPOINT ["python", "run.py"]