FROM python:3.9.18-bookworm
RUN apt-get update && apt-get upgrade -y
RUN yes | pip install --upgrade pip
RUN yes | pip install uropa==4.0.3 --exists-action i
COPY run_peak_annotation.py /run.py
ENTRYPOINT ["python", "run.py"]