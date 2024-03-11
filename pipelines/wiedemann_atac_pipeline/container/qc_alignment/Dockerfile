FROM python:3.13.0a4-bullseye
RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y openjdk-11-jre unzip r-base
RUN wget --no-verbose 'https://bitbucket.org/kokonech/qualimap/downloads/qualimap_v2.3.zip' && unzip 'qualimap_v2.3.zip' && rm 'qualimap_v2.3.zip'
ENV PATH="/qualimap_v2.3:${PATH}"
COPY qc_alignment.py /run.py
ENTRYPOINT ["python", "run.py"]