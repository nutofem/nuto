add_unit_test(MeshFem
    mechanics/mesh/MeshFemDofConvert
    mechanics/elements/ElementShapeFunctions)
target_link_libraries(MeshFem ${ANN_LIBRARIES})

add_unit_test(UnitMeshFem
    mechanics/mesh/MeshFem
    mechanics/elements/ElementShapeFunctions)
