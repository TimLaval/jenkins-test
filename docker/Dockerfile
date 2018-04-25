FROM ubuntu:16.04  

USER root

ARG DEBIAN_FRONTEND=noninteractive

ENV VTK_VERSION v6.3.0
ENV VTK_GIT https://gitlab.kitware.com/vtk/vtk.git

ENV N_CPUS 2

#Set update
RUN apt-get update && \
    apt-get -y install apt-utils locales && \
    echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen && \
    locale-gen en_US.utf8 && \
    /usr/sbin/update-locale LANG=en_US.UTF-8


# Set environment
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US.UTF-8
ENV TERM xterm

ENV HOME /work
ENV SOFT $HOME/soft
ENV BASHRC $HOME/.bashrc

# Create a non-priviledge user that will run the services
ENV BASICUSER basicuser
ENV BASICUSER_UID 1000

RUN useradd -m -d $HOME -s /bin/bash -N -u $BASICUSER_UID $BASICUSER && \
    mkdir $SOFT && \
    mkdir $HOME/.scripts && \
    mkdir $HOME/.nipype
USER $BASICUSER
WORKDIR $HOME


USER root
#Installation of tools
RUN apt-get -qq update
RUN apt-get -qq -y install \
cmake \
doxygen \
doxygen-gui \
g++ \
gcc \
gdb \
git \
gitk \
graphviz \
libgtest-dev \
meld \
sshpass && \
rm -rf /var/lib/apt/lists/* && \
update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-5   90 --slave /usr/bin/g++ g++ /usr/bin/g++-5 && \
update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.9 80 --slave /usr/bin/g++ g++ /usr/bin/g++-4.9 && \
apt-get build-dep vtk6

WORKDIR $SOFT
RUN \
    mkdir vtk && \
    cd vtk && \
    git clone $VTK_GIT -b $VTK_VERSION VTK && \
    mkdir build && \
    cd build && \
    cmake ../VTK && \
    make -j $N_CPUS && \
make install


#RUN cmake /usr/src/gtest/CMakeLists.txt

USER basicuser

RUN mkdir /repo
RUN mkdir /compil
RUN cd repo
# Clone the conf files into the docker container
# Have to add a user with his pass
#RUN git clone https://github.com/barco-healthcare/dermscan-ipi.git
RUN git clone https://github.com/TimLaval/gtest-jenkins.git


