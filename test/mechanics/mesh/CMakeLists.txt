add_unit_test(MeshFem
    mechanics/mesh/GeometryMeshFem
    mechanics/interpolation/InterpolationTrussLinear.cpp
    mechanics/interpolation/InterpolationTriangleLinear.cpp
    mechanics/interpolation/InterpolationQuadLinear.cpp
    mechanics/interpolation/InterpolationTriangleQuadratic.cpp
    mechanics/mesh/MeshFemDofConvert)

add_unit_test(GeometryMeshFem
    mechanics/interpolation/InterpolationTriangleLinear.cpp
    mechanics/interpolation/InterpolationTrussLinear.cpp
    mechanics/interpolation/InterpolationQuadLinear.cpp
    mechanics/interpolation/InterpolationBrickLinear.cpp
    )

add_unit_test(UnitMeshFem
    mechanics/mesh/GeometryMeshFem
    mechanics/interpolation/InterpolationTriangleLinear.cpp
    mechanics/interpolation/InterpolationTrussLinear.cpp
    mechanics/interpolation/InterpolationQuadLinear.cpp
    mechanics/interpolation/InterpolationBrickLinear.cpp
    )

create_symlink("meshes")

add_unit_test(MeshGmsh
    mechanics/mesh/GeometryMeshFem
    math/Legendre
    mechanics/interpolation/InterpolationTrussLinear.cpp
    mechanics/interpolation/InterpolationTrussLobatto.cpp
    mechanics/interpolation/InterpolationTrussQuadratic.cpp
    mechanics/interpolation/InterpolationTriangleLinear.cpp
    mechanics/interpolation/InterpolationTriangleQuadratic.cpp
    mechanics/interpolation/InterpolationQuadLinear.cpp
    mechanics/interpolation/InterpolationQuadLobatto.cpp
    mechanics/interpolation/InterpolationQuadSerendipity.cpp
    mechanics/interpolation/InterpolationTetrahedronLinear.cpp
    mechanics/interpolation/InterpolationTetrahedronQuadratic.cpp
    mechanics/interpolation/InterpolationBrickLinear.cpp
    mechanics/interpolation/InterpolationBrickLobatto.cpp
    mechanics/interpolation/InterpolationPrismLinear.cpp
    mechanics/interpolation/InterpolationPrismQuadratic.cpp
    mechanics/interpolation/InterpolationPyramidLinear.cpp)
