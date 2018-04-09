FROM debian

USER root

RUN apt-get update \
    && apt-get upgrade \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*
