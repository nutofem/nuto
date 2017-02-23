![alt text](doc/images/NuTo_logo.png "NuTo logo")

[![Build Status](https://travis-ci.org/nutofem/nuto.svg?branch=master)](https://travis-ci.org/nutofem/nuto)
[![codecov](https://codecov.io/gh/nutofem/nuto/branch/master/graph/badge.svg)](https://codecov.io/gh/nutofem/nuto)
[![Documentation Status](https://readthedocs.org/projects/nuto/badge/?version=master)](http://nuto.readthedocs.io/en/master/?badge=master)

What is it?
===========
A finite element library.

How do I install it?
====================
At your own peril.

First, you need some external dependencies

    sudo apt-get install git swig3.0 cmake doxygen python3-dev python3-numpy python3-instant\
    libboost-all-dev libeigen3-dev libopenblas-dev libmetis-dev libmumps-seq-dev libann-dev \
    libarpack2-dev libomp-dev gmsh

Then, you need to check out the source code

    git clone https://github.com/nutofem/nuto.git

Create a build directory and switch to it

    mkdir nuto/build && cd nuto/build

Run cmake

    cmake ..

Once this ran without errors, you can issue make (`-j4` for parallel building)

    make -j4

If you want to use the python module, and run all the tests, you need to add 
the module path to your environment

    export PYTHONPATH=<path/to/nuto/build>:$PYTHONPATH

In the end, run the test suite to see if all went well

    make test

Requirements
============

Name     | Version | Description
---------|---------|----------------
CMake    | ≥ 2.8   | Build system
Boost    | ≥ 1.54  | General purpose C++ libraries
Eigen    | ≥ 3.2   | C++ template library for linear algebra
Python   | ≥ 3.4   | Python programming language
NumPy    | ≥ 1.8   | Linear algebra in Python
OpenBLAS | ≥ 0.2.8 | Optimized BLAS library
METIS    | ≥ 5.1   | Serial graph partitioning and fill-reducing matrix ordering
MUMPS    |   4.10  | Parallel sparse direct solver
OpenMP   | ≥ 3.7   | Management of multiple OpenMP threads
ARPACK   | ≥ 3.1.5 | Solver for large scale eigenvalue problems
ANN      | ≥ 1.1.2 | Approximate nearest neighbour searching

Examples
========

![alt text](doc/images/crack_phase_field.png "Crack phase-field for a single edge notched tension test")

Can I see some code?
====================

Sure, take a look at [our documentation](https://nuto.readthedocs.io/en/master/).
