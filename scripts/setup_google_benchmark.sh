#!/bin/bash
git submodule init
git submodule update
cd external/benchmark
git checkout v1.3.0
cd ../..
