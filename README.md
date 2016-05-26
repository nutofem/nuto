NuTo
====

What is it?
===========
A finite element library.

How do I install it?
====================
At your own peril.

Can I see some code?
====================

Sure: @ref ipythonexample

As long as the new README is not ready, you can continue reading the old one here.
Most of it is very out of date (use git instead of svn, Mersenne is not needed
anymore, Eigen3 instead of 2, etc.), so tread carefully:

Old README
==========

NuTo is a project originally developed at the Institute of Structural Mechanics at the Bauhaus-University
in Weimar, Germany. It is OpenSource and the license is attached in LICENSE.txt

Instruction for compilation

1. checkout NuTo using 

		svn co https://nuto.svn.sourceforge.net/svnroot/nuto nuto

   or update using

   		svn up

2. checkout externals individually or from yesterday
  
        svn co svn+ssh://yesterday.bauing.uni-weimar.de/yesterday/slang/nuto/externals
		
Metis
-------

    cd metis
    tar xzf metis-4.0.tar.gz
    cd metis-4.0
    patch -p1 -i ../metis-4.0-cmake.patch
    patch -p1 -i ../metis-4.0-gcc44.patch (only required for gcc >= 4.4)
    mkdir objdir
    cd objdir
    cmake -DCMAKE_INSTALL_PREFIX=/opt/nuto ..
    make
    make install (eventually as root)

Mumps
-------

    cd mumps
    tar xvfz MUMPS_4.8.4.tar.gz
    cd MUMPS_4.8.4
    
    patch -p1 < ../MUMPS_4.8.4-cmake.patch
    mkdir objdir
    cd objdir
    cmake -DCMAKE_INSTALL_PREFIX=/opt/nuto ..
			
    make
    make install (eventually as root)
			
Mersenne Twister
----------------

    cd mersenne
    cd dSFMT-src-2.1
    patch -p1 < ../dSFMT-src-2.1-cmake.patch
    cmake -DCMAKE_INSTALL_PREFIX=/opt/nuto ..
    make
    make install (eventually as root)

Eigen2
-------

    cd eigen
    tar xzf 2.0.10.tar.gz
    cd eigen
    mkdir objdir
    cd objdir
    cmake -DCMAKE_INSTALL_PREFIX=/opt/nuto ..
    make
    make install (eventually as root)

BOOST
-------

    cd boost
    tar xzf boost_1_41_0.tar.gz
    cd boost_1_41_0
    ./bootstrap.sh --prefix=/opt/nuto
    ./bjam --layout=system  --with-serialization
    ./bjam --layout=system  --with-serialization --prefix=/opt/nuto install (eventually as root)
		
3. call cmake to build the Makefiles 

   out-of-source build in build directory 

        cmake mysourcedirectory -DNUTO_EXTERNAL_LIBRARIES_ROOT=/opt/nuto -DBOOST_ROOT=/opt/nuto

   or call cmake to generate an eclipse project

        cmake -G"Eclipse CDT4 - Unix Makefiles" -D CMAKE_BUILD_TYPE=Debug  -DNUTO_EXTERNAL_LIBRARIES_ROOT=/opt/nuto -DBOOST_ROOT=/opt/nuto ../nuto

4. call the build process 

		make
		
5. add module path to your environment

		export PYTHONPATH=/home/unger3/develop/nuto_build

6. execute test files

        make test

   a single test file (e.g. Brick8N) file can be executed

        ctest -V -R Brick8N
		
Requirements
------------

- BLAS
- LAPACK
- BOOST >= 1.36.0
- SWIG >=1.3.40 
- CMAKE >= 2.6.0
- Python >= 2.6
- Eigen >= 2.0.6
- Doxygen >= 1.5.8

Set Fortran-compiler with the envireonment variable FC or the cmake variable CMAKE_Fortran_COMPILER.

Compiling a debug version of python
-----------------------------------

    ./configure --without-pymalloc --enable-shared

uncomment in `Objects/obmalloc.c`:

    #define Py_USING_MEMORY_DEBUGGER

build

    make; make altinstall

Debugging an extension
----------------------

set breakpoint for gdb with

    br _PyImport_LoadDynamicModule

continue (`c`) until the corresponding module is loaded

    continue

update table of source files

    share

set breakpoint in your source file
continue debugging
