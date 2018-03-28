add_unit_test(DataArray)

add_unit_test(UnstructuredGrid
    visualize/Cell.cpp
    visualize/DataArray.cpp
    visualize/Point.cpp
    visualize/XMLWriter.cpp)

add_unit_test(VoronoiGeometries
    math/Legendre.cpp
    math/Quadrature.cpp
    mechanics/integrationtypes/IntegrationTypeTriangle.cpp
    mechanics/integrationtypes/IntegrationTypeTensorProduct.cpp
    mechanics/integrationtypes/IntegrationTypeTetrahedron.cpp)
