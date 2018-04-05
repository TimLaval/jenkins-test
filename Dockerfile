FROM debian:jessie

USER root

RUN useradd -d "/var/jenkins_home" -u 1000 -m -s /bin/bash jenkins
RUN mkdir -p /var/log/jenkins
RUN chown -R jenkins:jenkins /var/log/jenkins

VOLUME ["/var/log/jenkins"]

USER jenkins
CMD ["echo", "Data container for Jenkins"]

