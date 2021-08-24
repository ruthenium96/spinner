# syntax=docker/dockerfile:1
FROM ubuntu:20.04
ENV TZ=Europe/Moscow
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt update && apt install -y \
  cmake \
  git \
  lsb-release \
  wget \
  software-properties-common \
  gfortran \
  libopenblas-dev \
  liblapack-dev \
  libarpack2-dev \
  libsuperlu-dev \
  libarmadillo-dev \
  libboost-all-dev \
  && rm -rf /var/lib/apt/lists/*

RUN bash -c "$(wget -O - https://apt.llvm.org/llvm.sh)"

ENTRYPOINT ["tail", "-f", "/dev/null"]