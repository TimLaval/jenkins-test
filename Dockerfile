FROM debian

USER root

RUN apt-get install -y \
    build-essential \
cmake \
doxygen \
doxygen-gui \
g++ \
gdb \
git \
gitk \
libgl1-mesa-dev \
libgtest-dev \
libltdl-dev \
libncurses5-dev \
libxt-dev \
meld \
sshpass \
texlive \
utils

RUN cd /usr/src/gtest \
    && cmake CMakeLists.txt \
    && make


RUN apt-get update \
    && apt-get upgrade \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*
