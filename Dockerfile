FROM node:7-alpine

RUN curl -fsSLO https://get.docker.com/builds/Linux/x86_64/docker-17.04.0-ce.tgz \
    && tar xvzf docker-17.04.0-ce.tgz \
    && mv docker/docker /usr/local/bin \
    && rm -r docker docker-17.04.0-ce.tgz

RUN service docker stop \
    && nohup docker daemon -H tcp://0.0.0.0:2375 -H unix:///var/run/docker.sock &
