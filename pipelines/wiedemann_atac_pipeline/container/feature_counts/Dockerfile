FROM python:3.9.18-bookworm
RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y build-essential

# Installs featureCounts.
RUN wget --no-verbose 'https://altushost-swe.dl.sourceforge.net/project/subread/subread-2.0.6/subread-2.0.6-source.tar.gz' && \
    tar -xvf subread-2.0.6-source.tar.gz && \
    rm subread-2.0.6-source.tar.gz
WORKDIR /subread-2.0.6-source/src
RUN make -f Makefile.Linux
ENV PATH="/subread-2.0.6-source/bin:${PATH}"

WORKDIR /
COPY run_feature_counts.py /run.py
ENTRYPOINT ["python", "run.py"]