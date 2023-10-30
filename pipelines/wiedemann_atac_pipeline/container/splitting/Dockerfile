FROM python:3.11.6-alpine3.18
RUN apk -U upgrade && apk add --no-cache autoconf automake make gcc musl-dev perl bash zlib-dev bzip2-dev xz-dev curl-dev openssl-dev ncurses-dev g++ zlib unzip
RUN wget 'https://github.com/samtools/samtools/archive/refs/tags/1.18.zip' -O samtools.zip && unzip 'samtools.zip' && rm 'samtools.zip'
RUN wget 'https://github.com/samtools/htslib/releases/download/1.18/htslib-1.18.tar.bz2' -O htslib.tar.bz2 && tar -xf 'htslib.tar.bz2' && rm 'htslib.tar.bz2'

WORKDIR /htslib-1.18
RUN autoreconf -i
RUN ./configure
RUN make
RUN make install
WORKDIR /samtools-1.18
RUN autoreconf -i
RUN ./configure
RUN make
RUN make install

WORKDIR /
COPY splitting.py /run.py
ENTRYPOINT ["python", "run.py"]