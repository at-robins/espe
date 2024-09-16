FROM python:3.9.18-bookworm
RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y unzip wget
COPY genome_info /genome_info
COPY run.py /run.py
ENTRYPOINT ["python", "run.py"]