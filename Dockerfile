FROM debian:jessie

USER root

ENV DEBIAN_FRONTEND "noninteractive"

RUN apt-get -qq update
RUN apt-get -qq -y install \
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
texlive 


RUN useradd jenkins --shell /bin/bash --create-home
USER jenkins
