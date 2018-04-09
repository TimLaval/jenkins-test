FROM debian

COPY plugins.txt /usr/share/jenkins/ref/plugins.txt
RUN /usr/local/bin/install-plugins.sh < /usr/share/jenkins/ref/plugins.txt

RUN apt-get update \
&& apt-get clean \
&& rm -rf /var/lib/apt/lists/*
