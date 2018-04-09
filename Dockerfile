FROM debian:jessie

RUN useradd jenkins --shell /bin/bash --create-home
USER jenkins
