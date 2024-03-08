FROM python:3.9.18-bookworm
RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y gzip git
COPY run_blacklist.py /run.py
ENTRYPOINT ["python", "run.py"]