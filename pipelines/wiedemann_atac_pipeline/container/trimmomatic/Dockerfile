FROM python:3.11.6-alpine3.18
RUN apk -U upgrade && apk add --no-cache openjdk17-jre unzip
RUN wget 'https://github.com/usadellab/Trimmomatic/files/5854859/Trimmomatic-0.39.zip' && unzip 'Trimmomatic-0.39.zip' && rm 'Trimmomatic-0.39.zip'
COPY run_trimming.py /run_trimming.py
ENTRYPOINT ["python", "run_trimming.py"]