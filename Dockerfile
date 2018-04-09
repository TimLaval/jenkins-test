FROM debian:jessie

USER root

ENV DEBIAN_FRONTEND "noninteractive"
RUN apt-get install -y build-essential

RUN useradd jenkins --shell /bin/bash --create-home
USER jenkins
