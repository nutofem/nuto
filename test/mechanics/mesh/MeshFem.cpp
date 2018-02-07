#include "BoostUnitTest.h"
#include "mechanics/mesh/MeshFem.h"
#include "mechanics/mesh/MeshFemDofConvert.h"
#include "mechanics/interpolation/InterpolationTrussLinear.h"
#include "mechanics/interpolation/InterpolationTriangleLinear.h"
#include "mechanics/interpolation/InterpolationQuadLinear.h"
#include "mechanics/interpolation/InterpolationTriangleQuadratic.h"

void SetStuff(NuTo::MeshFem& m)
{
    m.Nodes[0].SetValue(0, 0);
}

NuTo::MeshFem DummyMesh(NuTo::DofType dofType)
{
    NuTo::MeshFem mesh;
    auto& interpolation = mesh.CreateInterpolation(NuTo::InterpolationTriangleLinear());

    auto& n0 = mesh.Nodes.Add(Eigen::Vector2d({1, 0}));
    auto& n1 = mesh.Nodes.Add(Eigen::Vector2d({2, 0}));
    auto& n2 = mesh.Nodes.Add(Eigen::Vector2d({0, 3}));

    auto& nd0 = mesh.Nodes.Add(42);
    auto& nd1 = mesh.Nodes.Add(4);
    auto& nd2 = mesh.Nodes.Add(6174);

    auto& e0 = mesh.Elements.Add({{{n0, n1, n2}, interpolation}});
    e0.AddDofElement(dofType, {{nd0, nd1, nd2}, interpolation});
    return mesh;
}

BOOST_AUTO_TEST_CASE(MeshAddStuff)
{
    NuTo::DofType d("Dof", 1);
    NuTo::MeshFem mesh = DummyMesh(d);

    auto& e0 = mesh.Elements[0];
    BoostUnitTest::CheckVector(e0.CoordinateElement().ExtractNodeValues(), std::vector<double>({1, 0, 2, 0, 0, 3}), 6);

    mesh.Nodes[0].SetValue(0, 4);
    BoostUnitTest::CheckVector(e0.CoordinateElement().ExtractNodeValues(), std::vector<double>({4, 0, 2, 0, 0, 3}), 6);

    NuTo::MeshFem meshMoved = std::move(mesh);
    meshMoved.Nodes[0].SetValue(0, 42);
    auto& e0FromMove = meshMoved.Elements[0];
    BoostUnitTest::CheckVector(e0FromMove.CoordinateElement().ExtractNodeValues(),
                               std::vector<double>({42, 0, 2, 0, 0, 3}), 6);
}

BOOST_AUTO_TEST_CASE(MeshNodeSelectionCoords)
{
    NuTo::DofType d("Dof", 1);
    NuTo::MeshFem mesh = DummyMesh(d);

    // selection of coordinate nodes
    {
        const auto& n = mesh.NodeAtCoordinate(Eigen::Vector2d(0, 3));
        BoostUnitTest::CheckEigenMatrix(n.GetValues(), Eigen::Vector2d(0, 3));
        BOOST_CHECK_THROW(mesh.NodeAtCoordinate(Eigen::Vector2d(0, 3.00001)), NuTo::Exception);
        BOOST_CHECK_NO_THROW(mesh.NodeAtCoordinate(Eigen::Vector2d(0, 3.00001), 1e-4));
    }


    // selection of dof nodes
    {
        const auto& n = mesh.NodeAtCoordinate(Eigen::Vector2d(0, 3), d);
        BOOST_CHECK_CLOSE(n.GetValues()[0], 6174, 1.e-10);
        BOOST_CHECK_THROW(mesh.NodeAtCoordinate(Eigen::Vector2d(0, 3.00001), d), NuTo::Exception);
        BOOST_CHECK_NO_THROW(mesh.NodeAtCoordinate(Eigen::Vector2d(0, 3.00001), d, 1e-4));
    }
}

BOOST_AUTO_TEST_CASE(MeshNodeSelectionAxis)
{
    NuTo::DofType d("Dof", 1);
    NuTo::MeshFem mesh = DummyMesh(d);

    {
        auto& nd0 = mesh.NodeAtCoordinate(Eigen::Vector2d(1, 0), d);
        auto& nd1 = mesh.NodeAtCoordinate(Eigen::Vector2d(2, 0), d);

        auto group0 = mesh.NodesAtAxis(NuTo::eDirection::Y, d);
        BOOST_CHECK_EQUAL(group0.Size(), 2);
        BOOST_CHECK(group0.Contains(nd0));
        BOOST_CHECK(group0.Contains(nd1));

        auto group1 = mesh.NodesAtAxis(NuTo::eDirection::Y, d, 2.);
        BOOST_CHECK(group1.Empty());
    }
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

BOOST_AUTO_TEST_CASE(MeshConvert)
{

    /* first create this:
     *
     * 2-----3
     * |\    |
     * | \   |
     * |  \  |
     * |   \ |
     * |    \|
     * 0-----1
     */
    NuTo::MeshFem mesh;
    auto& n0 = mesh.Nodes.Add(Eigen::Vector2d(0, 0));
    auto& n1 = mesh.Nodes.Add(Eigen::Vector2d(1, 0));
    auto& n2 = mesh.Nodes.Add(Eigen::Vector2d(0, 1));
    auto& n3 = mesh.Nodes.Add(Eigen::Vector2d(1, 1));

    auto& interpolation = mesh.CreateInterpolation(NuTo::InterpolationTriangleLinear());
    mesh.Elements.Add({{{n0, n1, n2}, interpolation}});
    mesh.Elements.Add({{{n1, n3, n2}, interpolation}});

    int expectedNumCoordinateNodes = 4;
    BOOST_CHECK_EQUAL(mesh.Nodes.Size(), expectedNumCoordinateNodes);


    // add linear dof type

    NuTo::DofType dof0("linear", 1);
    NuTo::AddDofInterpolation(&mesh, dof0, interpolation);

    int expectedNumDof0Nodes = expectedNumCoordinateNodes; // same interpolation

    BOOST_CHECK_EQUAL(mesh.Nodes.Size(), expectedNumCoordinateNodes + expectedNumDof0Nodes);
    BOOST_CHECK_NO_THROW(mesh.NodeAtCoordinate(Eigen::Vector2d(0, 0), dof0));


    /* Now add a quadratic interpolation
     *
     * 2--8--6
     * |\    |
     * | \   |
     * 5  4  7
     * |   \ |
     * |    \|
     * 0--3--1
     *
     * The numbering is not correct, but the total number of points is.
     */
    NuTo::DofType dof1("quadratic", 1);
    const auto& interpolationQuadratic = mesh.CreateInterpolation(NuTo::InterpolationTriangleQuadratic());
    NuTo::AddDofInterpolation(&mesh, dof1, interpolationQuadratic);

    int expectedNumDof1Nodes = 9;
    BOOST_CHECK_EQUAL(mesh.Nodes.Size(), expectedNumCoordinateNodes + expectedNumDof0Nodes + expectedNumDof1Nodes);
}

BOOST_AUTO_TEST_CASE(MeshConvertFromCoordinates)
{

    /* first create this:
     *
     * 2-----3-----5
     * |\    |     |
     * | \   |     |
     * |  \  |     |
     * |   \ |     |
     * |    \|     |
     * 0-----1-----4
     */
    NuTo::MeshFem mesh;
    auto& n0 = mesh.Nodes.Add(Eigen::Vector2d(0, 0));
    auto& n1 = mesh.Nodes.Add(Eigen::Vector2d(1, 0));
    auto& n2 = mesh.Nodes.Add(Eigen::Vector2d(0, 1));
    auto& n3 = mesh.Nodes.Add(Eigen::Vector2d(1, 1));
    auto& n4 = mesh.Nodes.Add(Eigen::Vector2d(2, 0));
    auto& n5 = mesh.Nodes.Add(Eigen::Vector2d(2, 1));

    auto& interpolationTriangle = mesh.CreateInterpolation(NuTo::InterpolationTriangleLinear());
    auto& interpolationQuad = mesh.CreateInterpolation(NuTo::InterpolationQuadLinear());
    mesh.Elements.Add({{{n0, n1, n2}, interpolationTriangle}});
    mesh.Elements.Add({{{n1, n4, n5, n3}, interpolationQuad}});

    int expectedNumCoordinateNodes = 6;
    BOOST_CHECK_EQUAL(mesh.Nodes.Size(), expectedNumCoordinateNodes);


    // add linear dof type

    NuTo::DofType dof0("isoparametric", 1);
    NuTo::AddDofInterpolation(&mesh, dof0);

    int expectedNumDof0Nodes = expectedNumCoordinateNodes; // same interpolation

    BOOST_CHECK_EQUAL(mesh.Nodes.Size(), expectedNumCoordinateNodes + expectedNumDof0Nodes);
    BOOST_CHECK_NO_THROW(mesh.NodeAtCoordinate(Eigen::Vector2d(0, 0), dof0));
}

BOOST_AUTO_TEST_CASE(MeshNodesTotalDof)
{
    NuTo::MeshFem mesh;
    auto& n0 = mesh.Nodes.Add(Eigen::Vector2d(0, 0));
    auto& n1 = mesh.Nodes.Add(Eigen::Vector2d(1, 0));
    auto& n2 = mesh.Nodes.Add(Eigen::Vector2d(0, 1));

    auto& interpolationTriangle = mesh.CreateInterpolation(NuTo::InterpolationTriangleLinear());
    auto& interpolationTruss = mesh.CreateInterpolation(NuTo::InterpolationTrussLinear());

    // Add coordinate elements
    auto& tri = mesh.Elements.Add({{{n0, n1, n2}, interpolationTriangle}});
    auto& line1 = mesh.Elements.Add({{{n0, n1}, interpolationTruss}});
    auto& line2 = mesh.Elements.Add({{{n1, n2}, interpolationTruss}});
    auto& line3 = mesh.Elements.Add({{{n2, n0}, interpolationTruss}});

    // Add a dof element
    NuTo::DofType dof1("dof1", 1);

    auto& nD0 = mesh.Nodes.Add(1.2);
    auto& nD1 = mesh.Nodes.Add(2.2);

    line1.AddDofElement(dof1, {{nD0, nD1}, interpolationTruss});

    BOOST_CHECK_EQUAL(mesh.NodesTotal(dof1).Size(), 2);
}
