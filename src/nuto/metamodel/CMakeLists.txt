SET(NuToMetamodel_SOURCES
  Metamodel.cpp
  MetamodelException.cpp
  MultipleLinearRegression.cpp
  NeuralNetwork.cpp
  TransferFunction.cpp
  Transformation.cpp
  MinMaxTransformation.cpp
  ZeroMeanUnitVarianceTransformation.cpp
  SupportPoints.cpp
)

SET(NuToMetamodel_SWIG_MODULES
  ModulMetaModel
)

# build static library
ADD_LIBRARY(NuToMetamodel STATIC ${NuToMetamodel_SOURCES})
IF(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
  SET_TARGET_PROPERTIES(NuToMetamodel PROPERTIES COMPILE_FLAGS "-fPIC")
ENDIF(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
# install static library
INSTALL(TARGETS NuToMetamodel ARCHIVE DESTINATION ${NUTO_LIBRARIES_INSTALL_PATH})

# build python modules
CONFIGURE_FILE(__init__.py ${CMAKE_CURRENT_BINARY_DIR}/__init__.py)
NUTO_INSTALL_PYTHON_FILE(__init__.py ${NUTO_PYTHON_MODULES_INSTALL_PATH}/nuto/metamodel)
FOREACH(NuToMetamodel_SWIG_MODULE ${NuToMetamodel_SWIG_MODULES})
  SET_SOURCE_FILES_PROPERTIES(${NuToMetamodel_SWIG_MODULE}.i PROPERTIES CPLUSPLUS ON)
  SET_SOURCE_FILES_PROPERTIES(${NuToMetamodel_SWIG_MODULE}.i PROPERTIES SWIG_FLAGS "${NuTo_SWIG_FLAGS}")
  SWIG_ADD_MODULE(${NuToMetamodel_SWIG_MODULE} python ${NuToMetamodel_SWIG_MODULE}.i)
  SWIG_LINK_LIBRARIES(${NuToMetamodel_SWIG_MODULE} NuToMetamodel NuToOptimize NuToMath NuToBase ${MERSENNE_LIBRARIES} ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${PYTHON_LIBRARIES})
  # check for unresolved symbols
  IF( NOT CYGWIN AND NOT WIN32 )
    SET_TARGET_PROPERTIES(${SWIG_MODULE_${NuToMetamodel_SWIG_MODULE}_REAL_NAME} PROPERTIES LINK_FLAGS -Wl,-z,defs)
  ENDIF( NOT CYGWIN AND NOT WIN32 )
  IF(MINGW)
    SET_TARGET_PROPERTIES(${SWIG_MODULE_${NuToMetamodel_SWIG_MODULE}_REAL_NAME} PROPERTIES LINK_FLAGS "-shared -Wl,--enable-auto-import")
  ENDIF(MINGW)
  NUTO_INSTALL_SWIG_PYTHON_MODULE(${NuToMetamodel_SWIG_MODULE} ${NUTO_PYTHON_MODULES_INSTALL_PATH}/nuto/metamodel)
ENDFOREACH(NuToMetamodel_SWIG_MODULE)
