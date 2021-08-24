# syntax=docker/dockerfile:1
FROM ubuntu:20.04
ENV TZ=Europe/Moscow
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt update

# Common dependencies
RUN apt install -y \
  lsb-release \
  software-properties-common \
  libboost-all-dev
  
# Aramadillo dependencies
RUN apt update && apt install -y \
  libopenblas-dev \
  liblapack-dev \
  libarpack2-dev \
  libsuperlu-dev \
  libarmadillo-dev

# wignerSymbols dependencies
RUN apt install -y \
    gfortran

# Tools
RUN apt install -y \
    cmake \
    git \
    wget

RUN bash -c "$(wget -O - https://apt.llvm.org/llvm.sh)"

RUN rm -rf /var/lib/apt/lists/*

ENTRYPOINT ["tail", "-f", "/dev/null"]