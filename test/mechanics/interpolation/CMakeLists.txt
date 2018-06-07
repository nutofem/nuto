add_unit_test(InterpolationTrussLinear)
add_unit_test(InterpolationTrussQuadratic)
add_unit_test(InterpolationTruss3rdOrder)
add_unit_test(InterpolationTruss4thOrder)

add_unit_test(InterpolationTetrahedronLinear
    mechanics/interpolation/InterpolationTrussLinear
    mechanics/interpolation/InterpolationTriangleLinear)

add_unit_test(InterpolationTetrahedronQuadratic
    mechanics/interpolation/InterpolationTrussQuadratic
    mechanics/interpolation/InterpolationTriangleQuadratic)

add_unit_test(InterpolationTriangleLinear
    mechanics/interpolation/InterpolationTrussLinear.cpp)
add_unit_test(InterpolationTriangleQuadratic
    mechanics/interpolation/InterpolationTrussQuadratic.cpp)
add_unit_test(InterpolationTriangle3rdOrder
    mechanics/interpolation/InterpolationTruss3rdOrder.cpp)
add_unit_test(InterpolationTriangle4thOrder
    mechanics/interpolation/InterpolationTruss4thOrder.cpp)

add_unit_test(InterpolationQuadLinear
    mechanics/interpolation/InterpolationTrussLinear.cpp)

add_unit_test(InterpolationQuadSerendipity
    mechanics/interpolation/InterpolationTrussQuadratic.cpp)

add_unit_test(InterpolationBrickLinear
    mechanics/interpolation/InterpolationTrussLinear.cpp
    mechanics/interpolation/InterpolationQuadLinear.cpp)

add_unit_test(InterpolationBrickQuadratic
    mechanics/interpolation/InterpolationTrussQuadratic.cpp
    mechanics/interpolation/InterpolationQuadSerendipity.cpp)

add_unit_test(InterpolationPrismLinear
    mechanics/interpolation/InterpolationTriangleLinear.cpp
    mechanics/interpolation/InterpolationQuadLinear.cpp
    mechanics/interpolation/InterpolationTrussLinear.cpp)

add_unit_test(InterpolationPrismQuadratic
    mechanics/interpolation/InterpolationTriangleQuadratic.cpp
    mechanics/interpolation/InterpolationTrussQuadratic.cpp
    mechanics/interpolation/InterpolationQuadSerendipity.cpp)

add_unit_test(InterpolationPyramidLinear
    mechanics/interpolation/InterpolationTriangleLinear.cpp
    mechanics/interpolation/InterpolationTrussLinear.cpp
    mechanics/interpolation/InterpolationQuadLinear.cpp)


add_unit_test(InterpolationTrussLobatto
    math/Legendre.cpp)
add_unit_test(InterpolationQuadLobatto
    mechanics/interpolation/InterpolationTrussLobatto.cpp
    math/Legendre.cpp)
add_unit_test(InterpolationBrickLobatto
    mechanics/interpolation/InterpolationTrussLobatto.cpp
    mechanics/interpolation/InterpolationQuadLobatto.cpp
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
