FROM python:3.9.18-bookworm
RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y samtools openjdk-17-jre r-base
RUN wget --no-verbose https://github.com/broadinstitute/picard/releases/download/3.1.1/picard.jar && chmod 774 /picard.jar
COPY duplicate_removal.py /run.py
ENTRYPOINT ["python", "run.py"]