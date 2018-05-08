FROM ubuntu

USER root

ARG DEBIAN_FRONTEND=noninteractive
# Host independent UID/GID mapping.
ENV USER=abcd USER_ID=1000 USER_GID=1000 WORK_DIR="/media" MKI_PATH="/usr/bin/mkisofs"
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

#Create the user which will run scripts.
RUN groupadd --gid ${USER_GID} ${USER}
RUN useradd -m -G sudo -s /bin/bash abcd --uid ${USER_ID} --gid ${USER_GID} --create-home
RUN echo "abcd ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers
RUN cat /etc/sudoers


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
lcov \
libgl1-mesa-dev \
libqt5x11extras5-dev \
libxt-dev \
meld \
ninja-build \
sshpass \
wget \
python3-dev \
python3-sip \
python3-sip-dev

RUN apt-get install -f

# Install CMake
ENV CMAKE_ARCHIVE_SHA256 33e4851d3219b720f4b64fcf617151168f1bffdf5afad25eb4b7f5f58cee3a08
ENV CMAKE_VERSION_MAJOR 3
ENV CMAKE_VERSION_MINOR 8
ENV CMAKE_VERSION_PATCH 2
ENV CMAKE_VERSION ${CMAKE_VERSION_MAJOR}.${CMAKE_VERSION_MINOR}.${CMAKE_VERSION_PATCH}
WORKDIR $SOFT
RUN wget https://cmake.org/files/v${CMAKE_VERSION_MAJOR}.${CMAKE_VERSION_MINOR}/cmake-${CMAKE_VERSION}-Linux-x86_64.tar.gz && \
  hash=$(sha256sum ./cmake-${CMAKE_VERSION}-Linux-x86_64.tar.gz | awk '{ print $1 }') && \
  [ $hash = "${CMAKE_ARCHIVE_SHA256}" ] && \
  tar -xzvf cmake-${CMAKE_VERSION}-Linux-x86_64.tar.gz && \
  rm cmake-${CMAKE_VERSION}-Linux-x86_64.tar.gz
ENV PATH=$SOFT/cmake-${CMAKE_VERSION}-Linux-x86_64/bin:${PATH}

#Install GoogleTest
RUN apt-get install -y libgtest-dev && \
    cd /usr/src/gtest && \
    cmake CMakeLists.txt && \
    make
 
# copy or symlink libgtest.a and libgtest_main.a to your /usr/lib folder
WORKDIR /usr/src/gtest
RUN cp *.a /usr/lib

#Install VTK
WORKDIR $SOFT
RUN mkdir VTK-build
RUN wget -nv -O- http://www.vtk.org/files/release/8.1/VTK-8.1.0.tar.gz | tar xz
WORKDIR $SOFT/VTK-8.1.0
    
RUN sed -i 's/\/\/#define\ GLX_GLXEXT_LEGACY/#define\ GLX_GLXEXT_LEGACY/g' Rendering/OpenGL/vtkXOpenGLRenderWindow.cxx

RUN apt-get update

RUN apt-get build-dep -y libvtk6-dev
WORKDIR $SOFT/VTK-build
RUN cmake -G Ninja \
  -DCMAKE_BUILD_TYPE:STRING=Release \
  -DBUILD_TESTING:BOOL=OFF \
  -DVTK_WRAP_PYTHON:BOOL=ON \
  -DVTK_WRAP_PYTHON_SIP:BOOL=ON \
  -DVTK_Group_Qt:BOOL=ON \
  -DVTK_QT_VERSION:STRING=5 \
  -DSIP_PYQT_DIR=/usr/share/sip/PyQt5/ \
  -DVTK_PYTHON_VERSION:STRING=3 \
  ../VTK-8.1.0
RUN ninja


WORKDIR $ABCD
RUN git clone https://github.com/TimLaval/gtest-jenkins.git
# Mkisofs support
COPY mkisofs /usr/bin
RUN chown abcd:abcd ${MKI_PATH}
RUN chmod u+x ${MKI_PATH}

COPY user_mapping.sh /
COPY runTest.sh /
COPY documentation.sh /
COPY coverage.sh /
RUN chmod u+x /user_mapping.sh
RUN chmod u+x /runTest.sh
RUN chmod u+x /documentation.sh
RUN chmod u+x /coverage.sh


USER ${USER}

