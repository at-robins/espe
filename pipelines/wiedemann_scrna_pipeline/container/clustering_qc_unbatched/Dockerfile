FROM python:3.9.18-bookworm
RUN apt-get update && apt-get upgrade -y

RUN yes | pip install --upgrade pip
RUN yes | pip install numpy==1.26.1 --exists-action i
RUN yes | pip install scanpy[leiden]==1.9.6 --exists-action i
RUN yes | pip install seaborn==0.12.2 --exists-action i

COPY run_clustering_unbatched.py /run.py
ENTRYPOINT ["python", "run.py"]