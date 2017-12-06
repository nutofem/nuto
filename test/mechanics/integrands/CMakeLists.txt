add_unit_test(NeumannBc
    mechanics/elements/ElementShapeFunctions.cpp)

add_unit_test(HistoryData
    math/Legendre.cpp
    mechanics/constraintsPde/Constraints.cpp
    mechanics/constraintsPde/ConstraintCompanion.cpp
    mechanics/dofs/DofNumbering.cpp
    mechanics/elements/ElementShapeFunctions.cpp
    mechanics/integrationtypes/IntegrationTypeBase.cpp
    mechanics/integrationtypes/IntegrationTypeTensorProduct.cpp
    mechanics/mesh/MeshFem.cpp
    mechanics/mesh/MeshFemDofConvert.cpp)

