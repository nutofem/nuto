set(MetamodelSources
    Metamodel.cpp
    NeuralNetwork.cpp
    TransferFunction.cpp
    Transformation.cpp
    MinMaxTransformation.cpp
    ZeroMeanUnitVarianceTransformation.cpp
    SupportPoints.cpp
    PolynomialLeastSquaresFitting.cpp
    )

sources_to_objects("${MetamodelSources}" MetamodelObjects)

add_library(Metamodel ${MetamodelObjects})
set_target_properties(Metamodel PROPERTIES
    OUTPUT_NAME NuToMetamodel
    )

target_link_libraries(Metamodel Optimize Math Base)

install(TARGETS Metamodel EXPORT NuToTargets
    LIBRARY DESTINATION lib
    )
