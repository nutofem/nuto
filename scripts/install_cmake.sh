#!/bin/bash
set -ev
# 3.5 is in Debian Stable and in 16.04, and has IMPORTED targets for boost
wget -q https://cmake.org/files/v3.5/cmake-3.5.2-Linux-x86_64.sh
sudo sh cmake-3.5.2-Linux-x86_64.sh -- --skip-license --prefix=/usr
