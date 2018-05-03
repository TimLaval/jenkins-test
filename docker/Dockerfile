FROM debian:jessie  

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
ENV ABCD $HOME/abcd
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
build-essential \ 
doxygen \
doxygen-gui \
g++ \
git \
gitk \
graphviz \
libgtest-dev \
libgl1-mesa-dev \
meld \
ninja-build \
sshpass \
wget

# Install CMake
ENV CMAKE_ARCHIVE_SHA256 33e4851d3219b720f4b64fcf617151168f1bffdf5afad25eb4b7f5f58cee3a08
ENV CMAKE_VERSION_MAJOR 3
ENV CMAKE_VERSION_MINOR 8
ENV CMAKE_VERSION_PATCH 2
ENV CMAKE_VERSION ${CMAKE_VERSION_MAJOR}.${CMAKE_VERSION_MINOR}.${CMAKE_VERSION_PATCH}
RUN cd $SOFT && \
  wget https://cmake.org/files/v${CMAKE_VERSION_MAJOR}.${CMAKE_VERSION_MINOR}/cmake-${CMAKE_VERSION}-Linux-x86_64.tar.gz && \
  hash=$(sha256sum ./cmake-${CMAKE_VERSION}-Linux-x86_64.tar.gz | awk '{ print $1 }') && \
  [ $hash = "${CMAKE_ARCHIVE_SHA256}" ] && \
  tar -xzvf cmake-${CMAKE_VERSION}-Linux-x86_64.tar.gz && \
  rm cmake-${CMAKE_VERSION}-Linux-x86_64.tar.gz
ENV PATH=$SOFT/cmake-${CMAKE_VERSION}-Linux-x86_64/bin:${PATH}

#Install VTK

RUN cd $SOFT && git clone --depth 1 git://vtk.org/VTK.git && \
    mkdir VTK-build && \
    cd VTK-build && \
    cmake \
        -G Ninja \
        -DCMAKE_INSTALL_PREFIX:PATH=/usr \
        -DBUILD_EXAMPLES:BOOL=OFF \
        -DBUILD_TESTING:BOOL=OFF \
        -DBUILD_SHARED_LIBS:BOOL=ON \
        -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=ON \
        -DITK_LEGACY_REMOVE:BOOL=ON \
        -DITK_BUILD_DEFAULT_MODULES:BOOL=OFF \
        -DITK_USE_SYSTEM_LIBRARIES:BOOL=OFF \
        -DModule_ITKCommon:BOOL=ON \
        -DModule_ITKIOXML:BOOL=ON \
        -DModule_ITKExpat:BOOL=ON \
        ../VTK && \
    ninja install && \
cd $SOFT && rm -rf VTK VTK-build
#RUN cmake /usr/src/gtest/CMakeLists.txt


RUN cd $HOME
RUN mkdir $ABCD
RUN mkdir $ABCD/repo
RUN mkdir $ABCD/compil
RUN cd $ABCD/repo
# Clone the conf files into the docker container
# Have to add a user with his pass
#RUN git clone https://github.com/barco-healthcare/dermscan-ipi.git
RUN git clone https://github.com/TimLaval/gtest-jenkins.git

USER basicuser

