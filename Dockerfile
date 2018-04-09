FROM debian:jessie

USER root
    
RUN useradd jenkins --shell /bin/bash --create-home
USER jenkins
