add_unit_test(InterpolationTrussLinear)

add_unit_test(InterpolationTetrahedronLinear)
add_unit_test(InterpolationTetrahedronQuadratic)

add_unit_test(InterpolationTriangleLinear)
add_unit_test(InterpolationTriangleQuadratic)

add_unit_test(InterpolationQuadLinear)
add_unit_test(InterpolationQuadSerendipity)

add_unit_test(InterpolationBrickLinear)

add_unit_test(InterpolationPrismLinear)
add_unit_test(InterpolationPrismQuadratic
    mechanics/interpolation/InterpolationTriangleQuadratic.cpp
    mechanics/interpolation/InterpolationTrussQuadratic.cpp)

add_unit_test(InterpolationPyramidLinear)


add_unit_test(InterpolationTrussLobatto
    mechanics/elements/SpectralShapeFunctions.cpp
    math/Legendre.cpp)
add_unit_test(InterpolationQuadLobatto
    mechanics/elements/SpectralShapeFunctions.cpp
    math/Legendre.cpp)
add_unit_test(InterpolationBrickLobatto
    mechanics/elements/SpectralShapeFunctions.cpp
    math/Legendre.cpp)

add_unit_test(InterpolationCompanion
    mechanics/elements/SpectralShapeFunctions.cpp
    math/Legendre.cpp
    )
