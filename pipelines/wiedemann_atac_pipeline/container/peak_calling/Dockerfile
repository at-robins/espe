FROM python:3.9.18-bookworm
RUN apt-get update && apt-get upgrade -y
RUN yes | pip install --upgrade pip
RUN yes | pip install macs3==3.0.1 --exists-action i
COPY run_peak_calling.py /run.py
ENTRYPOINT ["python", "run.py"]