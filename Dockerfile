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


RUN useradd jenkins --shell /bin/bash --create-home
RUN mkdir -p /data touch /data/x
RUN chown -R jenkins /data
VOLUME /data


RUN cmake /usr/src/gtest/CMakeLists.txt

USER jenkins

# Clone the conf files into the docker container
# Have to add a user with his pass
RUN git clone https://github.com/barco-healthcare/dermscan-ipi.git
