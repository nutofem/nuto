set(NuToMetamodel_SOURCES
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

# build static library
add_library(NuToMetamodel ${NuToMetamodel_SOURCES})
target_link_libraries(NuToMetamodel NuToOptimize NuToMath NuToBase)
if(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
    set_target_properties(NuToMetamodel PROPERTIES COMPILE_FLAGS "-fPIC")
endif()

# install library
install(TARGETS NuToMetamodel
    ARCHIVE DESTINATION ${NUTO_LIBRARIES_INSTALL_PATH}
    LIBRARY DESTINATION ${NUTO_LIBRARIES_INSTALL_PATH}
    )

# build python modules
configure_file(__init__.py ${CMAKE_CURRENT_BINARY_DIR}/__init__.py)
nuto_install_python_file(__init__.py ${NUTO_PYTHON_MODULES_INSTALL_PATH}/nuto/metamodel)
nuto_swig_module(ModulMetaModel ModulMetaModel.i nuto/metamodel NuToMetamodel)
