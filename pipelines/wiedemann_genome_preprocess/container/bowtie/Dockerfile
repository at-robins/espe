FROM python:3.11.6-alpine3.18
RUN apk -U upgrade && apk add --no-cache autoconf automake make gcc musl-dev perl bash zlib-dev bzip2-dev xz-dev curl-dev openssl-dev ncurses-dev g++ zlib unzip
RUN wget 'https://github.com/BenLangmead/bowtie2/archive/refs/tags/v2.5.1.zip' -O bowtie2.zip && unzip 'bowtie2.zip' && rm 'bowtie2.zip'

WORKDIR /bowtie2-2.5.1
RUN make
WORKDIR /
RUN mv /bowtie2-2.5.1 /bowtie2

COPY run.py /run.py
ENTRYPOINT ["python", "run.py"]