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
ENV TEMP $HOME/tempbuild
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
git \
gitk \
graphviz \
libgtest-dev \
meld \
sshpass \
wget

RUN mkdir /tmpbuild
WORKDIR /tmpbuild

# Download & build Tcl
# https://www.tcl.tk/doc/howto/compile.html#unix
RUN wget http://prdownloads.sourceforge.net/tcl/tcl8.6.6-src.tar.gz && tar -zxvf tcl8.6.6-src.tar.gz
RUN cd tcl8.6.6/unix && ./configure && make && make install

WORKDIR $TEMP
# Download & build Tk
# https://www.tcl.tk/doc/howto/compile.html
RUN wget http://prdownloads.sourceforge.net/tcl/tk8.6.6-src.tar.gz && tar -zxvf tk8.6.6-src.tar.gz
RUN cd tk8.6.6/unix && ./configure && make && make install

WORKDIR $TEMP
# Download & extract VTK
RUN wget http://www.vtk.org/files/release/7.0/VTK-7.0.0.tar.gz && tar -zxvf VTK-7.0.0.tar.gz

RUN mkdir /vtk-build2
RUN cd /vtk-build2/ && cmake \
  -DCMAKE_BUILD_TYPE:STRING=Release \
  -DBUILD_TESTING:BOOL=OFF \
  -DVTK_WRAP_PYTHON:BOOL=ON \
  -DVTK_WRAP_PYTHON_SIP:BOOL=ON \
  -DVTK_WRAP_TCL:BOOL=ON \
  -DVTK_PYTHON_VERSION:STRING=3 \
  -DVTK_USE_TK:BOOL=ON \
  /tmpbuild/VTK-7.0.0

# Build VTK
RUN cd /vtk-build2/ && make

#RUN cmake /usr/src/gtest/CMakeLists.txt

USER basicuser

RUN cd $HOME
RUN mkdir /repo
RUN mkdir /compil
RUN cd repo
# Clone the conf files into the docker container
# Have to add a user with his pass
#RUN git clone https://github.com/barco-healthcare/dermscan-ipi.git
RUN git clone https://github.com/TimLaval/gtest-jenkins.git


