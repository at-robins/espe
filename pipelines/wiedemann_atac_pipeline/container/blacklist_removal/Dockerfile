FROM python:3.9.18-bookworm
RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y bedtools samtools
COPY run_blacklist_removal.py /run.py
ENTRYPOINT ["python", "run.py"]