#!/bin/bash
set -ev
wget https://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O miniconda.sh
chmod +x miniconda.sh
bash miniconda.sh -b -p $HOME/miniconda
conda install --yes numpy
