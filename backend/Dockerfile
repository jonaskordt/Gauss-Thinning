FROM ubuntu:20.04

ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"
RUN apt-get update

RUN DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get -y install tzdata
RUN apt-get install -y wget
RUN apt-get install -y make build-essential libeigen3-dev

RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh 
RUN conda --version

RUN conda install -c conda-forge pybind11 igl websockets meshplot

RUN mkdir /code
WORKDIR /code

COPY third_party ./third_party
COPY src ./src

RUN ls third_party
RUN (cd third_party/DevelopableApproximationViaGaussImageThinning && make -j8)

WORKDIR src
ENTRYPOINT ["python", "main.py"]
