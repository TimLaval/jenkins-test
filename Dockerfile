FROM debian:jessie

USER root

ARG DEBIAN_FRONTEND "noninteractive"

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
gdb \
git \
gitk \
graphviz \
libgtest-dev \
meld \
sshpass

# Added some template
ADD doxygen/ doxygen/

# Go to working directory
#WORKDIR /doxygen/

WORKDIR $SOFT
RUN \
    mkdir vtk && \
    cd vtk && \
    git clone $VTK_GIT -b $VTK_VERSION VTK && \
    mkdir build && \
    cd build && \
    cmake -DCMAKE_BUILD_TYPE=Release \
          -DPYTHON_EXECUTABLE=/usr/bin/python2.7 \
          -DPYTHON_INCLUDE_DIR=/usr/include/python2.7 \
          -DPYTHON_INCLUDE_DIR=/usr/include/x86_64-linux-gnu/python2.7m \
          -DPYTHON_LIBRARY=/usr/lib/x86_64-linux-gnu/libpython2.7m.so.1 \
          ../VTK && \
    make -j $N_CPUS && \
make install


#RUN cmake /usr/src/gtest/CMakeLists.txt

USER basicuser

# Clone the conf files into the docker container
# Have to add a user with his pass
RUN git clone https://github.com/barco-healthcare/dermscan-ipi.git
