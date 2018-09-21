#include "BoostUnitTest.h"
#include "nuto/mechanics/mesh/GeometryMeshFem.h"
#include "nuto/mechanics/interpolation/InterpolationTrussLinear.h"
#include "nuto/mechanics/interpolation/InterpolationTriangleLinear.h"
#include "nuto/mechanics/interpolation/InterpolationQuadLinear.h"
#include "nuto/mechanics/interpolation/InterpolationTriangleQuadratic.h"

NuTo::GeometryMeshFem DummyMesh()
{
    NuTo::GeometryMeshFem mesh;
    auto& interpolation = mesh.CreateInterpolation(NuTo::InterpolationTriangleLinear());

    auto& n0 = mesh.AddNode(Eigen::Vector2d({1, 0}));
    auto& n1 = mesh.AddNode(Eigen::Vector2d({2, 0}));
    auto& n2 = mesh.AddNode(Eigen::Vector2d({0, 3}));

    mesh.AddElement({n0, n1, n2}, interpolation);

    return mesh;
}

BOOST_AUTO_TEST_CASE(MeshAddStuff)
{
    NuTo::GeometryMeshFem mesh = DummyMesh();

    auto& e0 = mesh.GetElement(0);
    BoostUnitTest::CheckVector(e0.ExtractCoordinates(), std::vector<double>({1, 0, 2, 0, 0, 3}), 6);

    mesh.ChangeNode(0, 0, 4);
    BoostUnitTest::CheckVector(e0.ExtractCoordinates(), std::vector<double>({4, 0, 2, 0, 0, 3}), 6);

    NuTo::GeometryMeshFem meshMoved = std::move(mesh);
    meshMoved.ChangeNode(0, 0, 42);
    auto& e0FromMove = meshMoved.GetElement(0);
    BoostUnitTest::CheckVector(e0FromMove.ExtractCoordinates(), std::vector<double>({42, 0, 2, 0, 0, 3}), 6);
}

BOOST_AUTO_TEST_CASE(MeshNodeSelectionCoords)
{
    NuTo::GeometryMeshFem mesh = DummyMesh();

    // selection of coordinate nodes
    {
        const auto& n = mesh.NodeAtCoordinate(Eigen::Vector2d(0, 3));
        BoostUnitTest::CheckEigenMatrix(n.GetCoordinates(), Eigen::Vector2d(0, 3));
        BOOST_CHECK_THROW(mesh.NodeAtCoordinate(Eigen::Vector2d(0, 3.00001)), NuTo::Exception);
        BOOST_CHECK_NO_THROW(mesh.NodeAtCoordinate(Eigen::Vector2d(0, 3.00001), 1e-4));
    }
}

BOOST_AUTO_TEST_CASE(MeshNodeSelectionAxis)
{
    NuTo::GeometryMeshFem mesh = DummyMesh();

    {
        auto& nd0 = mesh.NodeAtCoordinate(Eigen::Vector2d(1, 0));
        auto& nd1 = mesh.NodeAtCoordinate(Eigen::Vector2d(2, 0));

        auto group0 = mesh.NodesAtAxis(NuTo::eDirection::Y);
        BOOST_CHECK_EQUAL(group0.Size(), 2);
        BOOST_CHECK(group0.Contains(nd0));
        BOOST_CHECK(group0.Contains(nd1));

        auto group1 = mesh.NodesAtAxis(NuTo::eDirection::Y, 2.);
        BOOST_CHECK(group1.Empty());
    }
}

BOOST_AUTO_TEST_CASE(MeshSelectTotal)
{
    NuTo::GeometryMeshFem mesh = DummyMesh();
    auto nodes = mesh.NodesTotal();
    BOOST_CHECK_EQUAL(nodes.Size(), 3);
    auto elements = mesh.ElementsTotal();
    BOOST_CHECK_EQUAL(elements.Size(), 1);
}
