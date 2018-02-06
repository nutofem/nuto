add_unit_test(Cell
    mechanics/integrationtypes/IntegrationTypeBase.cpp
    mechanics/integrationtypes/IntegrationTypeEnum.cpp
    mechanics/elements/ElementShapeFunctions.cpp)

add_unit_test(Jacobian
    mechanics/elements/ElementShapeFunctions.cpp)

add_unit_test(CellData)

add_unit_test(SimpleAssembler)

add_unit_test(Matrix)

add_unit_test(CellIpData
    mechanics/elements/ElementShapeFunctions.cpp)
