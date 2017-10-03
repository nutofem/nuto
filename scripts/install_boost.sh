#!/bin/bash
set -ev
wget -q https://downloads.sourceforge.net/project/boost/boost/1.62.0/boost_1_62_0.tar.bz2
tar --bzip2 -xf boost_1_62_0.tar.bz2
cd boost_1_62_0
./bootstrap.sh --prefix=/usr --with-libraries=system,filesystem,test,mpi
echo "using mpi ;" >> project-config.jam
sudo ./b2 install -d0
