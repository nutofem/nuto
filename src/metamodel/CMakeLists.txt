set(NuToMetamodel_SOURCES
    Metamodel.cpp
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

# install library
install(TARGETS NuToMetamodel
    ARCHIVE DESTINATION ${NUTO_LIBRARIES_INSTALL_PATH}
    LIBRARY DESTINATION ${NUTO_LIBRARIES_INSTALL_PATH}
    )
