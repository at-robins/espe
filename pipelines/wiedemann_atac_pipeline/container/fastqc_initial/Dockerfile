FROM python:3.9.18-bookworm
RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y openjdk-17-jre unzip perl fontconfig fonts-dejavu
RUN wget 'https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip' && unzip 'fastqc_v0.12.1.zip' && rm 'fastqc_v0.12.1.zip'
COPY run_fastqc.py /run_fastqc.py
ENTRYPOINT ["python", "run_fastqc.py"]