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

    const auto& n = mesh.NodeAtCoordinate(Eigen::Vector2d(0, 3), d);
    BOOST_CHECK_CLOSE(n.GetValues()[0], 6174, 1.e-10);
    BOOST_CHECK_THROW(mesh.NodeAtCoordinate(Eigen::Vector2d(0, 0), d), NuTo::Exception);
}

BOOST_AUTO_TEST_CASE(MeshNodeSelectionAxis)
{
    NuTo::DofType d("Dof", 1);
    NuTo::MeshFem mesh = DummyMesh(d);
    
    auto& nd0 = mesh.NodeAtCoordinate(Eigen::Vector2d(1, 0), d);
    auto& nd1 = mesh.NodeAtCoordinate(Eigen::Vector2d(2, 0), d);

    auto group0 = mesh.NodeAtAxis(NuTo::eDirection::Y, 0., d);
    BOOST_CHECK_EQUAL(group0.Size(), 2);
    BOOST_CHECK(group0.Contains(nd0));
    BOOST_CHECK(group0.Contains(nd1));

    auto group1 = mesh.NodeAtAxis(NuTo::eDirection::Y, 2., d);
    BOOST_CHECK(group1.Empty());
}
