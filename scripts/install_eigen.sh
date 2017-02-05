#!/bin/bash
set -ev
wget http://bitbucket.org/eigen/eigen/get/3.3.2.tar.bz2
mkdir eigen && tar xf 3.3.2.tar.bz2 -C eigen --strip-components 1
mkdir eigen/build && cd eigen/build
cmake -DCMAKE_INSTALL_PREFIX=/usr ..
sudo make install
