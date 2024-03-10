FROM python:3.9.18-bookworm
RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y samtools bowtie2
COPY run_alignment.py /run_alignment.py
ENTRYPOINT ["python", "run_alignment.py"]