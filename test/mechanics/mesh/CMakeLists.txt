add_unit_test(MeshFem
    mechanics/mesh/MeshFemDofConvert
    mechanics/elements/ElementShapeFunctions)

add_unit_test(UnitMeshFem
    mechanics/mesh/MeshFem
    mechanics/elements/ElementShapeFunctions)

foreach(inputFile
    quad.msh
    quadNoncontiguous.msh
    )
create_symlink("${inputFile}")
endforeach()

add_unit_test(MeshGmsh
    mechanics/mesh/MeshFem
    mechanics/elements/ElementShapeFunctions)
