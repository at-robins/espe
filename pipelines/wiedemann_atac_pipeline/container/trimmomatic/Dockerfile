FROM python:3.9.18-bookworm
RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y openjdk-17-jre unzip
RUN wget 'https://github.com/usadellab/Trimmomatic/files/5854859/Trimmomatic-0.39.zip' && unzip 'Trimmomatic-0.39.zip' && rm 'Trimmomatic-0.39.zip'
COPY run_trimming.py /run_trimming.py
ENTRYPOINT ["python", "run_trimming.py"]