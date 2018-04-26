add_unit_test(NeumannBc)

add_unit_test(GradientDamage
    math/Legendre.cpp
    math/Quadrature.cpp
    mechanics/elements/SpectralShapeFunctions.cpp)

add_unit_test(HistoryData
    math/Legendre.cpp
    math/Quadrature.cpp
    mechanics/integrationtypes/IntegrationTypeTensorProduct.cpp)

