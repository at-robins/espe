FROM python:3.9.18-bookworm
RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y samtools

RUN yes | pip install --upgrade pip
RUN yes | pip install pysam==0.22.0 --exists-action i

COPY splitting.py /run.py
ENTRYPOINT ["python", "run.py"]