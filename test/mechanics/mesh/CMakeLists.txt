add_unit_test(MeshFem
    mechanics/mesh/MeshFemDofConvert
    mechanics/elements/ElementShapeFunctions)

add_unit_test(UnitMeshFem
    mechanics/mesh/MeshFem
    mechanics/elements/ElementShapeFunctions)

foreach(inputFile
    binary.msh
    quad.msh
    quadNoncontiguous.msh
    quad1.msh
    quadSerendipity2.msh
    triangle1.msh
    triangle2.msh
    line1.msh
    tet1.msh
    hex1.msh
    prism1.msh
    pyramid1.msh
    )
create_symlink("${inputFile}")
endforeach()

add_unit_test(MeshGmsh
    mechanics/mesh/MeshFem
    mechanics/elements/ElementShapeFunctions)
