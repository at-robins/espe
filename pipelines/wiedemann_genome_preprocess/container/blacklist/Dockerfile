FROM python:3.11.6-alpine3.18
RUN apk -U upgrade && apk add --no-cache gzip git
COPY run.py /run.py
ENTRYPOINT ["python", "run.py"]