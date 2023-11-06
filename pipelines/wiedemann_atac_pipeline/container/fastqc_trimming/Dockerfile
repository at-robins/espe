FROM python:3.11.6-alpine3.18
RUN apk -U upgrade && apk add --no-cache openjdk17-jre unzip perl fontconfig ttf-dejavu
RUN wget 'https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip' && unzip 'fastqc_v0.12.1.zip' && rm 'fastqc_v0.12.1.zip'
COPY run_fastqc.py /run_fastqc.py
ENTRYPOINT ["python", "run_fastqc.py"]