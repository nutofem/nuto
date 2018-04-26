add_unit_test(MeshFem
    mechanics/mesh/MeshFemDofConvert)

add_unit_test(UnitMeshFem
    mechanics/mesh/MeshFem)

foreach(inputFile
    binary.msh
    brickL2.msh
    quad.msh
    quadL2.msh
    quadNoncontiguous.msh
    quad1.msh
    quadSerendipity2.msh
    triangle1.msh
    triangle2.msh
    line1.msh
    line2.msh
    tet1.msh
    hex1.msh
    prism1.msh
    pyramid1.msh
    )
create_symlink("${inputFile}")
endforeach()

add_unit_test(MeshGmsh
    mechanics/mesh/MeshFem
    math/Legendre
    mechanics/elements/SpectralShapeFunctions.cpp)
