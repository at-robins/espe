FROM python:3.9.18-bookworm
RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y autoconf automake make gcc g++ pkg-config libhdf5-dev llvm-14-dev libssl-dev libcurl4-openssl-dev less locales wget ca-certificates fonts-texgyre libharfbuzz-dev libfribidi-dev r-base gfortran libglpk-dev libgmp-dev libxml2-dev

RUN yes | pip install --upgrade pip
RUN yes | pip install numpy==1.26.1 --exists-action i
RUN yes | pip install scanpy[leiden]==1.9.6 --exists-action i
RUN yes | pip install seaborn==0.12.2 --exists-action i
RUN yes | pip install pathvalidate==3.2.0 --exists-action i

# Installs Rust.
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- --default-toolchain=1.76.0 -y
ENV PATH="/root/.cargo/bin:${PATH}"

# Installs the cluster stability analysis tool.
RUN mkdir /tmp_downloads && \
    git -C /tmp_downloads clone https://github.com/at-robins/leiden-optimisation.git && \
    git -C /tmp_downloads/leiden-optimisation checkout 73e0e35
WORKDIR /tmp_downloads/leiden-optimisation
RUN cargo build --release && \
    cp /tmp_downloads/leiden-optimisation/target/release/leiden-optimisation /usr/bin && \
    rm -r /tmp_downloads

WORKDIR /
COPY run_cluster_relation_tree.py /run.py
ENTRYPOINT ["python", "run.py"]