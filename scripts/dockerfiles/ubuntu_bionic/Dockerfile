# 18.04
FROM ubuntu:bionic

# generic C++ stuff, we somehow need ssh to make MPI work (needed by gmsh)
RUN apt-get update && apt-get install -y --no-install-recommends g++ clang cmake make libboost-filesystem-dev libboost-mpi-dev libboost-test-dev lcov curl sudo ssh ca-certificates && rm -rf /var/lib/apt/lists/*
# scientific stuff
RUN apt-get update && apt-get install -y --no-install-recommends libeigen3-dev libiomp-dev python3-dev python3-numpy libopenblas-dev libmetis-dev libmumps-seq-dev libarpack2-dev gmsh && rm -rf /var/lib/apt/lists/*

#fix missing /dev/fd from https://github.com/jbbarth/docker-ruby/commit/1916309122b7c04be4c01c46910471fc1d8176c6
RUN test -e /dev/fd || ln -s /proc/self/fd /dev/fd

# setup non root user named nuto with sudo rights, following the instructions from
# https://stackoverflow.com/questions/25845538/using-sudo-inside-a-docker-container
RUN useradd nuto && echo "nuto:nuto" | chpasswd && adduser nuto sudo && mkdir -p /home/nuto && chown -R nuto:nuto /home/nuto
USER nuto

# create the source and build directories
RUN mkdir /home/nuto/source && mkdir /home/nuto/build

WORKDIR /home/nuto/build
