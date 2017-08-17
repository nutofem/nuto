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

create_nuto_module(Metamodel "${MetamodelSources}")
target_link_libraries(Metamodel Optimize Math Base)
