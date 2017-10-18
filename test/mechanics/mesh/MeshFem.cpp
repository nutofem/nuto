#include "BoostUnitTest.h"
#include "mechanics/mesh/MeshFem.h"
#include "mechanics/interpolation/InterpolationTriangleLinear.h"

NuTo::MeshFem DummyMesh(NuTo::DofType dofType)
{
    NuTo::MeshFem mesh;
    auto& interpolationCoords = mesh.CreateInterpolation(NuTo::InterpolationTriangleLinear(2));
    auto& interpolationDof = mesh.CreateInterpolation(NuTo::InterpolationTriangleLinear(dofType.GetNum()));

    auto& n0 = mesh.Nodes.Add(Eigen::Vector2d({1, 0}));
    auto& n1 = mesh.Nodes.Add(Eigen::Vector2d({2, 0}));
    auto& n2 = mesh.Nodes.Add(Eigen::Vector2d({0, 3}));

    auto& nd0 = mesh.Nodes.Add(42);
    auto& nd1 = mesh.Nodes.Add(4);
    auto& nd2 = mesh.Nodes.Add(6174);

    auto& e0 = mesh.Elements.Add({{{n0, n1, n2}, interpolationCoords}});
    e0.AddDofElement(dofType, {{nd0, nd1, nd2}, interpolationDof});
    return mesh;
}

BOOST_AUTO_TEST_CASE(MeshAddStuff)
{
    NuTo::DofType d("Dof", 1);
    NuTo::MeshFem mesh = DummyMesh(d);

    auto& e0 = mesh.Elements.front();
    BoostUnitTest::CheckVector(e0.CoordinateElement().ExtractNodeValues(), std::vector<double>({1, 0, 2, 0, 0, 3}), 6);
}

BOOST_AUTO_TEST_CASE(MeshNodeSelectionCoords)
{
    NuTo::DofType d("Dof", 1);
    NuTo::MeshFem mesh = DummyMesh(d);

    // selection of coordinate nodes
    {
        const auto& n = mesh.NodeAtCoordinate(Eigen::Vector2d(0, 3));
        BoostUnitTest::CheckEigenMatrix(n.GetValues(), Eigen::Vector2d(0, 3));
        BOOST_CHECK_THROW(mesh.NodeAtCoordinate(Eigen::Vector2d(0, 0)), NuTo::Exception);
    }


    // selection of dof nodes
    {
        const auto& n = mesh.NodeAtCoordinate(Eigen::Vector2d(0, 3), d);
        BOOST_CHECK_CLOSE(n.GetValues()[0], 6174, 1.e-10);
        BOOST_CHECK_THROW(mesh.NodeAtCoordinate(Eigen::Vector2d(0, 0), d), NuTo::Exception);
    }
}

BOOST_AUTO_TEST_CASE(MeshNodeSelectionAxis)
{
    NuTo::DofType d("Dof", 1);
    NuTo::MeshFem mesh = DummyMesh(d);

    auto& nd0 = mesh.NodeAtCoordinate(Eigen::Vector2d(1, 0), d);
    auto& nd1 = mesh.NodeAtCoordinate(Eigen::Vector2d(2, 0), d);

    auto group0 = mesh.NodesAtAxis(NuTo::eDirection::Y, d);
    BOOST_CHECK_EQUAL(group0.Size(), 2);
    BOOST_CHECK(group0.Contains(nd0));
    BOOST_CHECK(group0.Contains(nd1));

    auto group1 = mesh.NodesAtAxis(NuTo::eDirection::Y, d, 2.);
    BOOST_CHECK(group1.Empty());
}


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

    NuTo::UnitMeshFem::Transform(&mesh, f);

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
