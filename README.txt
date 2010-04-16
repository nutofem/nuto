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
			
dex@vnet01:~/Uni/nuto/externals/metis/metis-4.0/objdir$ make
Scanning dependencies of target metis
[  1%] Building C object Lib/CMakeFiles/metis.dir/balance.c.o
In file included from /home/dex/Uni/nuto/externals/metis/metis-4.0/Lib/metis.h:36,
                 from /home/dex/Uni/nuto/externals/metis/metis-4.0/Lib/balance.c:16:
/home/dex/Uni/nuto/externals/metis/metis-4.0/Lib/proto.h:462: error: conflicting types for ‘__log2’
/usr/include/bits/mathcalls.h:145: note: previous declaration of ‘__log2’ was here
make[2]: *** [Lib/CMakeFiles/metis.dir/balance.c.o] Error 1
make[1]: *** [Lib/CMakeFiles/metis.dir/all] Error 2
make: *** [all] Error 2

			
			
		 Mumps
		-------
		    cd mumps
		    tar xvfz MUMPS_4.8.4.tar.gz
			cd MUMPS_4.8.4
			
			patch -p1 < ../MUMPS_4.8.4-cmake.patch
			mkdir objdir
			cd objdir
			cmake -DCMAKE_INSTALL_PREFIX=/opt/nuto ..
			
braucht metis
	
	dex@vnet01:~/Uni/nuto/externals/mumps/MUMPS_4.8.4/objdir$ cmake -DCMAKE_INSTALL_PREFIX=/opt/nuto ..
-- The Fortran compiler identification is GNU
-- Check for working Fortran compiler: /usr/bin/f95
-- Check for working Fortran compiler: /usr/bin/f95 -- works
-- Checking whether /usr/bin/f95 supports Fortran 90
-- Checking whether /usr/bin/f95 supports Fortran 90 -- yes
-- Searching for Metis library ...
CMake Error at /usr/share/cmake-2.6/Modules/FindPackageHandleStandardArgs.cmake:57 (MESSAGE):
  Could NOT find Metis (missing: Metis_LIBRARIES Metis_INCLUDE_DIR)
Call Stack (most recent call first):
  CMakeModules/FindMetis.cmake:30 (find_package_handle_standard_args)
  CMakeLists.txt:58 (FIND_PACKAGE)

			make
			make install (eventually as root)
			
		 Mersenne Twister
		-------
		    cd mersenne
FEHLT: 	tar xzvf dSFMT-src-2.1.tar.gz
			cd dSFMT-src-2.1
            patch -p1 < ../dSFMT-src-2.1-cmake.patch
			cmake -DCMAKE_INSTALL_PREFIX=/opt/nuto ..
MUSSTE cMakeLists.txt eins hoeher kopieren
			make
			make install (eventually as root)
MUSSTE dSMFT.h eins hoeher kopieren		



dex@vnet01:~/Uni/nuto/externals/mersenne/dSFMT-src-2.1$ make
Scanning dependencies of target dSFMT
[100%] Building C object CMakeFiles/dSFMT.dir/dSFMT.c.o
In file included from /home/dex/Uni/nuto/externals/mersenne/dSFMT-src-2.1/dSFMT-params.h:4,
                 from /home/dex/Uni/nuto/externals/mersenne/dSFMT-src-2.1/dSFMT.c:17:
/home/dex/Uni/nuto/externals/mersenne/dSFMT-src-2.1/dSFMT.h:39:4: warning: #warning "DSFMT_MEXP is not defined. I assume DSFMT_MEXP is 19937."
Linking C shared library libdSFMT.so
[100%] Built target dSFMT


dex@vnet01:~/Uni/nuto/externals/mersenne/dSFMT-src-2.1$ sudo make install
[100%] Built target dSFMT
Install the project...
-- Install configuration: ""
CMake Error at cmake_install.cmake:36 (FILE):
  file INSTALL cannot find file
  "/home/dex/Uni/nuto/externals/mersenne/dSFMT.h" to install.


make: *** [install] Error 1

		
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
LIEF NICHT
BESSER:  build.sh		

		
		 BOOST
		-------
            cd boost
            tar xzf boost_1_41_0.tar.gz
            cd boost_1_41_0
            ./bootstrap.sh --prefix=/opt/nuto
            ./bjam --layout=system  --with-serialization
            ./bjam --layout=system  --with-serialization --prefix=/opt/nuto install (eventually as root)
		
		only for developers : to create a patch  
		    diff -rupN original_directory new_directory > patchfile, where new_directory now adds the CMakeFile
		



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
  BLAS
  LAPACK
  BOOST >= 1.36.0
  SWIG >=1.3.40 
  CMAKE >= 2.6.0
  Python >= 2.6
  Eigen >= 2.0.6
  Doxygen >= 1.5.8

Set Fortran-compiler with the envireonment variable FC or the cmake variable CMAKE_Fortran_COMPILER.

Compiling a debug version of python
   ./configure --without-pymalloc --enable-shared
   uncomment in Objects/obmalloc.c #define Py_USING_MEMORY_DEBUGGER
   make; make altinstall

Debugging an extension
   set breakpoint for gdb with
      br _PyImport_LoadDynamicModule
   continue (c) until the corresponding module is loaded
      continue
   update table of source files
       share
   set breakpoint in your source file
   continue debugging
