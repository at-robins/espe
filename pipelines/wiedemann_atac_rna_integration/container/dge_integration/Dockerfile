FROM python:3.9.18-bookworm
RUN apt-get update && apt-get upgrade -y

RUN yes | pip install --upgrade pip
RUN yes | pip install numpy==1.26.4 --exists-action i
RUN yes | pip install pandas==2.2.2 --exists-action i
RUN yes | pip install seaborn==0.13.2 --exists-action i
RUN yes | pip install scipy==1.13.1 --exists-action i
RUN yes | pip install fastcluster==1.2.6 --exists-action i

COPY run_dge_integration.py /run.py
ENTRYPOINT ["python", "/run.py"]