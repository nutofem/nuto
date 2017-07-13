#!/bin/bash
set -ev
# 2.8.5 is the version in debian stable
wget http://gmsh.info/bin/Linux/gmsh-2.8.5-Linux64.tgz
tar xf gmsh-2.8.5-Linux64.tgz
sudo cp gmsh-2.8.5-Linux/bin/gmsh /usr/bin
