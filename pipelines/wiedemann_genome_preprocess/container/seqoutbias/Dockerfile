FROM python:3.9.18-bookworm
RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y samtools genometools build-essential

# Installs Rust.
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- --default-toolchain=1.73.0 -y
ENV PATH="/root/.cargo/bin:${PATH}"

# Installs wigToBigWig
WORKDIR /
RUN git clone git://genome-source.soe.ucsc.edu/kent.git
WORKDIR /kent/src
RUN make utils
RUN find /root/bin -name "wigToBigWig" | xargs cp -t /usr/bin
WORKDIR /
RUN rm -r /kent

RUN wigToBigWig
COPY run_seqout.py /run.py
ENTRYPOINT ["python", "run.py"]