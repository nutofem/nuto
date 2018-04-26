add_unit_test(Cell
    mechanics/integrationtypes/IntegrationTypeEnum.cpp
    mechanics/interpolation/InterpolationQuadLinear.cpp
    )

add_unit_test(Jacobian
    mechanics/interpolation/InterpolationTrussLinear.cpp
    mechanics/interpolation/InterpolationTrussQuadratic.cpp
    mechanics/interpolation/InterpolationTriangleLinear.cpp
    mechanics/interpolation/InterpolationQuadLinear.cpp
    mechanics/interpolation/InterpolationTetrahedronLinear.cpp
    )

add_unit_test(CellData)

add_unit_test(DifferentialOperators)

add_unit_test(SimpleAssembler)

add_unit_test(Matrix)

add_unit_test(CellIpData)
