FROM python:3.9.18-bookworm
RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y bedtools

# Installs Rust.
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- --default-toolchain=1.92.0 -y
ENV PATH="/root/.cargo/bin:${PATH}"

# Installs the consensus peak tool.
RUN mkdir /tmp_downloads && \
    git -C /tmp_downloads clone https://github.com/at-robins/Gipfelkreuzer.git && \
    git -C /tmp_downloads/Gipfelkreuzer checkout f195087
WORKDIR /tmp_downloads/Gipfelkreuzer
RUN cargo build --release && \
    cp /tmp_downloads/Gipfelkreuzer/target/release/Gipfelkreuzer /usr/bin && \
    rm -r /tmp_downloads
WORKDIR /

COPY run_peak_merging.py /run.py
ENTRYPOINT ["python", "run.py"]