FROM debian

USER root

RUN cd /usr/src/gtest \
    && cmake CMakeLists.txt \
    && make


RUN apt-get update \
    && apt-get upgrade \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*
