add_unit_test(NeumannBc
    mechanics/interpolation/InterpolationTrussLinear.cpp
    mechanics/interpolation/InterpolationTriangleLinear.cpp
    mechanics/interpolation/InterpolationTriangleQuadratic.cpp
    )

add_unit_test(GradientDamage
    math/Legendre.cpp
    math/Quadrature.cpp)

add_unit_test(HistoryData
    math/Legendre.cpp
    math/Quadrature.cpp
    mechanics/interpolation/InterpolationTrussLinear.cpp
    mechanics/integrationtypes/IntegrationTypeTensorProduct.cpp)

