FROM lhofmann/paraview-vtk-ubuntu16:latest

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


#Installation of tools

RUN apt-get -qq update
RUN apt-get -qq upgrade
RUN apt-get -qq -y install \
cppcheck \
curl \
doxygen \
doxygen-gui \
g++ \
graphviz \
lcov \
libjpeg62 \
libnetcdf11 \
libgl1-mesa-dev \
libqt5x11extras5-dev \
libxt-dev \
libtiff5-dev \
libavcodec-dev \
libavformat-dev \
libswscale-dev \
libv4l-dev \
libxvidcore-dev \
libx264-dev \
libgtk-3-dev \
libatlas-base-dev \
gfortran \
meld \
sshpass \
wget \
pkg-config \
python3 \
python-opencv \
software-properties-common

RUN apt-get install -f

#Install GoogleTest
RUN apt-get install -y libgtest-dev && \
    cd /usr/src/gtest && \
    cmake CMakeLists.txt && \
    make
 
# copy or symlink libgtest.a and libgtest_main.a to your /usr/lib folder
WORKDIR /usr/src/gtest
RUN cp *.a /usr/lib


WORKDIR $SOFT
RUN apt-get install -qq python3-tk
RUN curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
RUN python3 get-pip.py
RUN pip3 install matplotlib
RUN pip3 install opencv-python
RUN pip3 install Pillow
 
WORKDIR $ABCD
# Mkisofs support
COPY mkisofs /usr/bin
RUN chown abcd:abcd ${MKI_PATH}
RUN chmod u+x ${MKI_PATH}

COPY user_mapping.sh /
COPY Doxyfile $ABCD

COPY analysis.sh $ABCD
COPY build_abcd.sh $ABCD
COPY documentation.sh $ABCD
COPY segmentation_tool.sh $ABCD
COPY testing_and_coverage.sh $ABCD

RUN chmod u+x /user_mapping.sh
#RUN chmod u+x $ABCD/exec_test_cov_doc.sh

