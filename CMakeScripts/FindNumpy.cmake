execute_process(
COMMAND /usr/bin/env python -c "import numpy; print numpy.get_include()"
OUTPUT_VARIABLE PYTHON_SYS_PATH
)
string(STRIP ${PYTHON_SYS_PATH} PYTHON_SYS_PATH)
FIND_PATH(NUMPY_INCLUDE_PATH arrayobject.h ${PYTHON_SYS_PATH}/numpy)

IF(NOT NUMPY_I_PATH AND NOT DEFINED $ENV{NUMPY_I_PATH})
  SET(NUMPY_SEARCH_DIRS  /usr/lib/python2.7/dist-packages/instant/swig)
  find_path(NUMPY_I_PATH "numpy.i" HINTS ${NUMPY_SEARCH_DIRS})
#  SET(NUMPY_I_PATH       /usr/lib/python2.7/dist-packages/instant/swig)  
ENDIF(NOT NUMPY_I_PATH AND NOT DEFINED $ENV{NUMPY_I_PATH})

set(rt "rt-NOTFOUND")
find_file(rt "numpy.i" ${NUMPY_I_PATH} NO_DEFAULT_PATH)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(NUMPY DEFAULT_MSG NUMPY_INCLUDE_PATH)
