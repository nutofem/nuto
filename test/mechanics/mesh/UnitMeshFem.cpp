#include "BoostUnitTest.h"
#include "mechanics/mesh/UnitMeshFem.h"

void Check2DMesh(NuTo::MeshFem mesh)
{
    BOOST_CHECK_EQUAL(mesh.Nodes.size(), 3 * 8);

    BOOST_CHECK_NO_THROW(mesh.NodeAtCoordinate(Eigen::Vector2d(0., 0.)));
    BOOST_CHECK_NO_THROW(mesh.NodeAtCoordinate(Eigen::Vector2d(1., 1.)));

    BOOST_CHECK_NO_THROW(mesh.NodeAtCoordinate(Eigen::Vector2d(1. / 2., 1. / 7.)));
    BOOST_CHECK_NO_THROW(mesh.NodeAtCoordinate(Eigen::Vector2d(1. / 2., 5. / 7.)));


    auto f = [](Eigen::VectorXd coords) // transforms mesh to (42,4) -- (44, 11)
    {
        Eigen::VectorXd newCoords(2);
        newCoords[0] = 42 + coords[0] * 2;
        newCoords[1] = 4 + coords[1] * 7;
        return newCoords;
    };

    NuTo::MeshFem transformedMesh = NuTo::UnitMeshFem::Transform(std::move(mesh), f);

    BOOST_CHECK_NO_THROW(mesh.NodeAtCoordinate(Eigen::Vector2d(42., 4.)));
    BOOST_CHECK_NO_THROW(mesh.NodeAtCoordinate(Eigen::Vector2d(44., 11.)));

    BOOST_CHECK_NO_THROW(mesh.NodeAtCoordinate(Eigen::Vector2d(43., 5.)));
    BOOST_CHECK_NO_THROW(mesh.NodeAtCoordinate(Eigen::Vector2d(43., 9.)));
}

BOOST_AUTO_TEST_CASE(MeshQuad)
{
    auto mesh = NuTo::UnitMeshFem::CreateQuads(2, 7);
    BOOST_CHECK_EQUAL(mesh.Elements.size(), 2 * 7);
    Check2DMesh(mesh);
}

BOOST_AUTO_TEST_CASE(MeshTriangle)
{
    auto mesh = NuTo::UnitMeshFem::CreateTriangles(2, 7);
    BOOST_CHECK_EQUAL(mesh.Elements.size(), 2 * 7 * 2);
    Check2DMesh(mesh);
}
