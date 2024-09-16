FROM python:3.9.18-bookworm
RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y samtools genometools build-essential

# Installs faCount.
RUN mkdir /tmp_downloads && git -C /tmp_downloads clone git://genome-source.soe.ucsc.edu/kent.git
RUN git -C /tmp_downloads/kent checkout 29a377ab11ee9503b5e1823ad44b18322967e921
WORKDIR /tmp_downloads/kent/src
RUN make utils
RUN find /root/bin -name "faCount" | xargs cp -t /usr/bin

# Clean-up.
WORKDIR /
RUN rm -r /tmp_downloads

COPY run_size.py /run.py
ENTRYPOINT ["python", "run.py"]