FROM debian:jessie

USER root
RUN apt-get update && apt-get install -y libelf1

RUN useradd jenkins --shell /bin/bash --create-home
USER jenkins
