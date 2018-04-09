FROM debian:jessie

USER root
RUN apt-get update \
    && apt-get install -y libelf1 \
    && apt-get install -y build-essential \
    && apt-get install -y cmake \
    && apt-get install -y doxygen \
    && apt-get install -y doxygen-gui \
    && apt-get install -y g++ \
    && apt-get install -y gdb \
    && apt-get install -y git \
    && apt-get install -y gitk \
    && apt-get install -y libgl1-mesa-dev \
    && apt-get install -y libgtest-dev \
    && apt-get install -y libltdl-dev \
    && apt-get install -y libncurses5-dev \
    && apt-get install -y libxt-dev \
    && apt-get install -y meld \
    && apt-get install -y sshpass \
    && apt-get install -y texlive \
    && apt-get install §y utils


RUN useradd jenkins --shell /bin/bash --create-home
USER jenkins
