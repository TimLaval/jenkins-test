FROM debian:jessie

USER root

ENV DEBIAN_FRONTEND "noninteractive"

RUN apt-get -qq update
RUN apt-get -qq -y install build-essential

RUN useradd jenkins --shell /bin/bash --create-home
USER jenkins
