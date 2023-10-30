FROM python:3.11.6-alpine3.18
RUN apk -U upgrade && apk add --no-cache curl make git cairo-dev pango-dev gcc g++
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- --default-toolchain=1.73.0 -y
RUN source $HOME/.cargo/env
RUN wget 'https://github.com/genometools/genometools/archive/refs/tags/v1.6.4.tar.gz' -O genometools.tar.gz && tar -xvf 'genometools.tar.gz' && rm 'genometools.tar.gz'


# RUN wget 'https://github.com/samtools/htslib/releases/download/1.18/htslib-1.18.tar.bz2' -O htslib.tar.bz2 && tar -xf 'htslib.tar.bz2' && rm 'htslib.tar.bz2'

WORKDIR /genometools-1.6.4
# RUN autoreconf -i
# RUN ./configure
RUN make -j4
RUN make install

WORKDIR /
RUN git clone git://genome-source.soe.ucsc.edu/kent.git
RUN apk add --no-cache mariadb-connector-c-dev bash shadow
WORKDIR /kent/src
# Sets bash as default shell to prevent errors during makefile script execution.
RUN sed -i 's@/bin/ash@bin/bash@g' /etc/passwd
# Skips building the hg directory since we only need the wigToBigWig utility.
# Also prevents a current error while building the hg utilities.
RUN sed -i 's@cd hg && ${MAKE} utils@echo "skipping hg"@g' makefile
RUN make

# WORKDIR /samtools-1.18
# RUN autoreconf -i
# RUN ./configure
# RUN make
# RUN make install

WORKDIR /
RUN wigToBigWig
COPY run.py /run.py
ENTRYPOINT ["python", "run.py"]