FROM python:3.11.6-alpine3.18
RUN apk -U upgrade && apk add --no-cache autoconf automake make gcc g++ pkgconfig hdf5 hdf5-dev llvm14 llvm14-dev llvm14-static llvm14-libs

# RUN yes | pip install numpy==1.26.1 scanpy==1.9.6 seaborn==0.12.2 --exists-action i
RUN yes | pip install --upgrade pip
RUN yes | pip install numpy==1.26.1 --exists-action i
RUN ls -alth /usr/bin/
# RUN yes | LLVM_CONFIG=/usr/bin/llvm14-config pip install llvmlite==0.41.1 --exists-action i
RUN yes | LLVM_CONFIG=/usr/bin/llvm14-config pip install 'scanpy[leiden]'==1.9.6 --exists-action i
RUN yes | pip install seaborn==0.12.2 --exists-action i

WORKDIR /
COPY run_preprocessing.py /run.py
ENTRYPOINT ["python", "run.py"]