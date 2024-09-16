FROM python:3.9.18-bookworm
RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y bowtie2
COPY run_bowtie.py /run.py
ENTRYPOINT ["python", "run.py"]