FROM debian:jessie

USER root

ENV DEBIAN_FRONTEND "noninteractive"

RUN apt-get -qq update
RUN apt-get -qq -y install \
cmake \
doxygen \
doxygen-gui \
g++ \
gdb \
git \
gitk \
libgtest-dev \
meld \
sshpass

RUN cmake /usr/src/gtest/CMakeLists.txt


RUN useradd jenkins --shell /bin/bash --create-home
USER jenkins
