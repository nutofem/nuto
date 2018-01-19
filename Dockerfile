FROM ubuntu:xenial

RUN apt-get update && apt-get dist-upgrade -y
# Toolchain?
RUN apt-get install -y g++ clang cmake
# boost
RUN apt-get install -y libboost-filesystem-dev libboost-system-dev libboost-mpi-dev libboost-test-dev

# Nuto deps
RUN apt-get install -y libeigen3-dev libiomp-dev python3-dev python3-numpy libopenblas-dev libmetis-dev libmumps-seq-dev libann-dev libarpack2-dev gmsh 

# coverage and documentation
RUN apt-get install -y lcov curl 

#fix missing /dev/fd from https://github.com/jbbarth/docker-ruby/commit/1916309122b7c04be4c01c46910471fc1d8176c6
RUN test -e /dev/fd || ln -s /proc/self/fd /dev/fd

# setup non root user named nuto with sudo rights, following the instructions from
# https://stackoverflow.com/questions/25845538/using-sudo-inside-a-docker-container
RUN apt-get install -y sudo
RUN useradd nuto && echo "nuto:nuto" | chpasswd && adduser nuto sudo
RUN mkdir -p /home/nuto && chown -R nuto:nuto /home/nuto
USER nuto

# create the source directory /home/nuto/source
RUN mkdir /home/nuto/source
RUN mkdir /home/nuto/build

# the build directory will be in /home/nuto/source/build
WORKDIR /home/nuto/build
