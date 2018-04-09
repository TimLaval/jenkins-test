FROM debian:jessie

USER root

ENV DEBIAN_FRONTEND "noninteractive"

RUN useradd jenkins --shell /bin/bash --create-home
USER jenkins
