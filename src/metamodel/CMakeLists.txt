set(NuToMetamodelSources
    Metamodel.cpp
    NeuralNetwork.cpp
    TransferFunction.cpp
    Transformation.cpp
    MinMaxTransformation.cpp
    ZeroMeanUnitVarianceTransformation.cpp
    SupportPoints.cpp
    PolynomialLeastSquaresFitting.cpp
    )

sources_to_objects("${NuToMetamodelSources}" NuToMetamodelObjects)

# build static library
add_library(NuToMetamodel ${NuToMetamodelObjects})
target_link_libraries(NuToMetamodel NuToOptimize NuToMath NuToBase)

install(TARGETS NuToMetamodel
    LIBRARY DESTINATION lib
    )
