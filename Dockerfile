FROM debian:jessie

USER root
RUN apt-get update \
    && apt-get upgrade
    
RUN useradd jenkins --shell /bin/bash --create-home
USER jenkins
