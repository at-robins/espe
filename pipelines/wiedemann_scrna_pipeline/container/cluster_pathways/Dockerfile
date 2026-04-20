FROM python:3.12.12-trixie
RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y autoconf automake make gcc g++ pkg-config libhdf5-dev llvm-19-dev libssl-dev libcurl4-openssl-dev less locales wget ca-certificates fonts-texgyre libharfbuzz-dev libfribidi-dev r-base gfortran libglpk-dev libgmp-dev libxml2-dev

RUN yes | pip install --upgrade pip
RUN yes | pip install numpy==2.4.2 --exists-action i
RUN yes | pip install scanpy[leiden]==1.12 --exists-action i
RUN yes | pip install seaborn==0.13.2 --exists-action i
RUN yes | pip install pathvalidate==3.3.1 --exists-action i
RUN yes | pip install liana==1.7.1 --exists-action i
RUN yes | pip install decoupler==2.1.6 --exists-action i

# Downloads the required pathway network files.
COPY download_pathways.py /download_pathways.py
RUN python /download_pathways.py

# Downloads the required GMT files.
RUN mkdir /gmts
RUN wget -O /gmts/reactome_human.gmt https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Hs/c2.cp.reactome.v2026.1.Hs.symbols.gmt
RUN wget -O /gmts/reactome_mouse.gmt https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Mm/m2.cp.reactome.v2026.1.Mm.symbols.gmt
RUN wget -O /gmts/go_human.gmt https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Hs/c5.all.v2026.1.Hs.symbols.gmt
RUN wget -O /gmts/go_mouse.gmt https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Mm/m5.all.v2026.1.Mm.symbols.gmt
RUN wget -O /gmts/hallmark_human.gmt https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Hs/h.all.v2026.1.Hs.symbols.gmt
RUN wget -O /gmts/hallmark_mouse.gmt https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Mm/mh.all.v2026.1.Mm.symbols.gmt

COPY run_cluster_pathways.py /run.py
ENTRYPOINT ["python", "/run.py"]