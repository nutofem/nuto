add_unit_test(InterpolationTrussLinear)
add_unit_test(InterpolationTrussQuadratic)
add_unit_test(InterpolationTruss3rdOrder)
add_unit_test(InterpolationTruss4thOrder)

add_unit_test(InterpolationTetrahedronLinear)
add_unit_test(InterpolationTetrahedronQuadratic)

add_unit_test(InterpolationTriangleLinear)
add_unit_test(InterpolationTriangleQuadratic)
add_unit_test(InterpolationTriangle3rdOrder)
add_unit_test(InterpolationTriangle4thOrder)

add_unit_test(InterpolationQuadLinear)
add_unit_test(InterpolationQuadSerendipity)

add_unit_test(InterpolationBrickLinear)
add_unit_test(InterpolationBrickQuadratic)

add_unit_test(InterpolationPrismLinear)
add_unit_test(InterpolationPrismQuadratic
    mechanics/interpolation/InterpolationTriangleQuadratic.cpp
    mechanics/interpolation/InterpolationTrussQuadratic.cpp)

add_unit_test(InterpolationPyramidLinear)


add_unit_test(InterpolationTrussLobatto
    math/Legendre.cpp)
add_unit_test(InterpolationQuadLobatto
    mechanics/interpolation/InterpolationTrussLobatto.cpp
    math/Legendre.cpp)
add_unit_test(InterpolationBrickLobatto
    mechanics/interpolation/InterpolationTrussLobatto.cpp
    math/Legendre.cpp)

add_unit_test(InterpolationCompanion
    mechanics/interpolation/InterpolationTrussLinear.cpp
    mechanics/interpolation/InterpolationTrussQuadratic.cpp
    mechanics/interpolation/InterpolationTrussLobatto.cpp
    mechanics/interpolation/InterpolationTriangleLinear.cpp
    mechanics/interpolation/InterpolationTriangleQuadratic.cpp
    mechanics/interpolation/InterpolationQuadLinear.cpp
    mechanics/interpolation/InterpolationQuadSerendipity.cpp
    mechanics/interpolation/InterpolationQuadLobatto.cpp
    mechanics/interpolation/InterpolationTetrahedronLinear.cpp
    mechanics/interpolation/InterpolationTetrahedronQuadratic.cpp
    mechanics/interpolation/InterpolationBrickLinear.cpp
    mechanics/interpolation/InterpolationBrickLobatto.cpp
    mechanics/interpolation/InterpolationPrismLinear.cpp
    mechanics/interpolation/InterpolationPrismQuadratic.cpp
    mechanics/interpolation/InterpolationPyramidLinear.cpp
    math/Legendre.cpp
    )
