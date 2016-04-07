SET(NuToMetamodel_SOURCES
  Metamodel.cpp
  MultipleLinearRegression.cpp
  NeuralNetwork.cpp
  TransferFunction.cpp
  Transformation.cpp
  MinMaxTransformation.cpp
  ZeroMeanUnitVarianceTransformation.cpp
  SupportPoints.cpp
  PolynomialLeastSquaresFitting.cpp
)

SET(NuToMetamodel_SWIG_MODULES
  ModulMetaModel
)

# build static library
ADD_LIBRARY(NuToMetamodel ${NuToMetamodel_SOURCES})
TARGET_LINK_LIBRARIES(NuToMetamodel NuToOptimize NuToMath NuToBase ${MERSENNE_LIBRARIES})
IF(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
  SET_TARGET_PROPERTIES(NuToMetamodel PROPERTIES COMPILE_FLAGS "-fPIC")
ENDIF(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
# install library
INSTALL(TARGETS NuToMetamodel
	ARCHIVE DESTINATION ${NUTO_LIBRARIES_INSTALL_PATH}
	LIBRARY DESTINATION ${NUTO_LIBRARIES_INSTALL_PATH})

# build python modules
CONFIGURE_FILE(__init__.py ${CMAKE_CURRENT_BINARY_DIR}/__init__.py)
NUTO_INSTALL_PYTHON_FILE(__init__.py ${NUTO_PYTHON_MODULES_INSTALL_PATH}/nuto/metamodel)
FOREACH(NuToMetamodel_SWIG_MODULE ${NuToMetamodel_SWIG_MODULES})
  NUTO_SWIG_MODULE(${NuToMetamodel_SWIG_MODULE} ${NuToMetamodel_SWIG_MODULE}.i nuto/metamodel NuToMetamodel)
ENDFOREACH(NuToMetamodel_SWIG_MODULE)
