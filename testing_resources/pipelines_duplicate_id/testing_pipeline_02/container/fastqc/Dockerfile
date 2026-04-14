# This is just a dummy docker file for testing.
FROM alpine:3.17.3
RUN apk -U upgrade && apk add --no-cache openjdk17-jre unzip perl fontconfig ttf-dejavu
RUN wget 'https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip' && unzip 'fastqc_v0.12.1.zip' && rm 'fastqc_v0.12.1.zip'
RUN chmod 755 /FastQC/fastqc
ENTRYPOINT ["perl" , "--", "/FastQC/fastqc", "--outdir=/output"]