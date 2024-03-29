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
  
# Armadillo dependencies
RUN apt install -y \
  libopenblas-dev \
  liblapack-dev \
  libarpack2-dev \
  libsuperlu-dev \
  libarmadillo-dev

# wignerSymbols dependencies
RUN apt install -y \
    gfortran

# OpenMP:
RUN apt install -y libomp-dev

# Tools
RUN apt install -y \
    cmake \
    git \
    wget \
    python3-pip \
    cppcheck

# Linter for cpp
RUN python3 -m pip install cpplint

# Install clang 13 compiler
RUN bash -c "$(wget -O - https://apt.llvm.org/llvm.sh)"

# Add clang-format
RUN apt install -y clang-format

RUN rm -rf /var/lib/apt/lists/*

ENTRYPOINT ["tail", "-f", "/dev/null"]