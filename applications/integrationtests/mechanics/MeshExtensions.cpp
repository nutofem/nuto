#include "BoostUnitTest.h"

#include "nuto/mechanics/mesh/MeshCompanion.h"
#include "nuto/mechanics/nodes/NodeSimple.h"
#include "nuto/mechanics/mesh/MeshFem.h"
#include "nuto/mechanics/mesh/UnitMeshFem.h"
#include "nuto/mechanics/interpolation/InterpolationTriangleQuadratic.h"

using namespace NuTo;

BOOST_AUTO_TEST_CASE(Triangle2ndOrderMeshAddEdges)
{
    MeshFem mesh;

    auto& nd0 = mesh.Nodes.Add(Eigen::Vector2d(0., 0.));
    auto& nd1 = mesh.Nodes.Add(Eigen::Vector2d(1., 0.));
    auto& nd2 = mesh.Nodes.Add(Eigen::Vector2d(1., 1.));
    auto& nd3 = mesh.Nodes.Add(Eigen::Vector2d(0., 1.));

    auto& ndA1 = mesh.Nodes.Add(Eigen::Vector2d(0.5, 0.));
    auto& ndA2 = mesh.Nodes.Add(Eigen::Vector2d(1.0, 0.5));
    auto& ndA3 = mesh.Nodes.Add(Eigen::Vector2d(0.5, 1.0));
    auto& ndA4 = mesh.Nodes.Add(Eigen::Vector2d(0., 0.5));
    auto& ndA5 = mesh.Nodes.Add(Eigen::Vector2d(0.5, 0.5));

    InterpolationTriangleQuadratic ipol;

    mesh.Elements.Add({{{nd0, nd1, nd2, ndA1, ndA2, ndA5}, ipol}});
    mesh.Elements.Add({{{nd2, nd3, nd0, ndA3, ndA4, ndA5}, ipol}});

    Group<ElementCollectionFem> allElements = mesh.ElementsTotal();

    auto edgesAll = AddEdgeElements(&mesh, allElements);

    BOOST_CHECK_EQUAL(edgesAll.Size(), 5);
}

BOOST_AUTO_TEST_CASE(BrickLinearMeshAddEdgesAndFaces)
{
    std::cout << "Create a test mesh" << std::endl;

    MeshFem mesh = UnitMeshFem::CreateBricks(2, 1, 1);

    Group<ElementCollectionFem> allElements = mesh.ElementsTotal();

    auto edgesAll = AddEdgeElements(&mesh, allElements);
    BOOST_CHECK_EQUAL(edgesAll.Size(), 20);

    auto facesAll = AddFaceElements(&mesh, allElements);
    BOOST_CHECK_EQUAL(facesAll.Size(), 11);

    auto facesOrientedAll = AddFaceElements(&mesh, allElements, true);
    BOOST_CHECK_EQUAL(facesOrientedAll.Size(), 12);
}
