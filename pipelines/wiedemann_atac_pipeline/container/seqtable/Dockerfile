FROM python:3.9.18-bookworm
RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y samtools genometools build-essential

# Installs Rust.
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- --default-toolchain=1.73.0 -y
ENV PATH="/root/.cargo/bin:${PATH}"

# Installs wigToBigWig, bigWigToBedGraph and bigWigMerge.
RUN mkdir /tmp_downloads && git -C /tmp_downloads clone git://genome-source.soe.ucsc.edu/kent.git
WORKDIR /tmp_downloads/kent/src
RUN make utils
RUN find /root/bin -name "wigToBigWig" | xargs cp -t /usr/bin && \
    find /root/bin -name "bigWigToBedGraph" | xargs cp -t /usr/bin && \
    find /root/bin -name "bigWigMerge" | xargs cp -t /usr/bin

# Installs seqOutBias.
RUN git -C /tmp_downloads clone https://github.com/guertinlab/seqOutBias.git && git -C /tmp_downloads/seqOutBias checkout 1446ebc
WORKDIR /tmp_downloads/seqOutBias
RUN cargo build --release && cp /tmp_downloads/seqOutBias/target/release/seqOutBias /usr/bin

# Clean-up.
WORKDIR /
RUN rm -r /tmp_downloads

COPY run_seqtable.py /run.py
ENTRYPOINT ["python", "run.py"]