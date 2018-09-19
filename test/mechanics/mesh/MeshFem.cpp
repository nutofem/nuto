#include "BoostUnitTest.h"
#include "nuto/mechanics/mesh/MeshFem.h"
#include "nuto/mechanics/mesh/GeometryMeshFem.h"
#include "nuto/mechanics/mesh/MeshFemDofConvert.h"
#include "nuto/mechanics/interpolation/InterpolationTrussLinear.h"
#include "nuto/mechanics/interpolation/InterpolationTriangleLinear.h"
#include "nuto/mechanics/interpolation/InterpolationQuadLinear.h"
#include "nuto/mechanics/interpolation/InterpolationTriangleQuadratic.h"

void SetStuff(NuTo::MeshFem& m)
{
    m.Nodes[0].SetValue(0, 0);
}

NuTo::MeshFem DummyMesh(NuTo::GeometryMeshFem& geoMesh, NuTo::DofType dofType)
{
    NuTo::MeshFem mesh(geoMesh);
    auto& interpolation = mesh.CreateInterpolation(NuTo::InterpolationTriangleLinear());

    auto& n0 = geoMesh.CoordinateNodes.Add(Eigen::Vector2d({1, 0}));
    auto& n1 = geoMesh.CoordinateNodes.Add(Eigen::Vector2d({2, 0}));
    auto& n2 = geoMesh.CoordinateNodes.Add(Eigen::Vector2d({0, 3}));

    auto& nd0 = mesh.Nodes.Add(42);
    auto& nd1 = mesh.Nodes.Add(4);
    auto& nd2 = mesh.Nodes.Add(6174);

    auto& cElm = geoMesh.Elements.Add({{n0, n1, n2}, interpolation});

    auto& e0 = mesh.Elements.Add(cElm);
    e0.AddDofElement(dofType, {{nd0, nd1, nd2}, interpolation});

    // Add another element without a interpolation for `dofType`. This must not
    // trigger exceptions in node select methods, when they try to find
    // the missing DofElement.

    auto& cElm2 = geoMesh.Elements.Add({{n0, n1, n2}, interpolation});
    mesh.Elements.Add(cElm2);

    return mesh;
}

BOOST_AUTO_TEST_CASE(AllocateInstances)
{
    NuTo::DofType d("Dof", 1);
    NuTo::GeometryMeshFem geoMesh;
    NuTo::MeshFem mesh = DummyMesh(geoMesh, d);
    mesh.AllocateDofInstances(d, 10);
    for (auto node : mesh.NodesTotal(d))
    {
        BOOST_CHECK_EQUAL(node.GetNumInstances(), 10);
    }
}

BOOST_AUTO_TEST_CASE(MeshNodeSelectionCoords)
{
    NuTo::DofType d("Dof", 1);
    NuTo::GeometryMeshFem geoMesh;
    NuTo::MeshFem mesh = DummyMesh(geoMesh, d);

    // selection of coordinate nodes
    {
        const auto& n = geoMesh.NodeAtCoordinate(Eigen::Vector2d(0, 3));
        BoostUnitTest::CheckEigenMatrix(n.GetCoordinates(), Eigen::Vector2d(0, 3));
        BOOST_CHECK_THROW(geoMesh.NodeAtCoordinate(Eigen::Vector2d(0, 3.00001)), NuTo::Exception);
        BOOST_CHECK_NO_THROW(geoMesh.NodeAtCoordinate(Eigen::Vector2d(0, 3.00001), 1e-4));
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
    NuTo::GeometryMeshFem geoMesh;
    NuTo::MeshFem mesh = DummyMesh(geoMesh, d);

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
        auto& nd0 = geoMesh.NodeAtCoordinate(Eigen::Vector2d(1, 0));
        auto& nd1 = geoMesh.NodeAtCoordinate(Eigen::Vector2d(2, 0));

        auto group0 = geoMesh.NodesAtAxis(NuTo::eDirection::Y);
        BOOST_CHECK_EQUAL(group0.Size(), 2);
        BOOST_CHECK(group0.Contains(nd0));
        BOOST_CHECK(group0.Contains(nd1));

        auto group1 = geoMesh.NodesAtAxis(NuTo::eDirection::Y, 2.);
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
    NuTo::GeometryMeshFem geoMesh;
    NuTo::MeshFem mesh(geoMesh);
    auto& n0 = geoMesh.CoordinateNodes.Add(Eigen::Vector2d(0, 0));
    auto& n1 = geoMesh.CoordinateNodes.Add(Eigen::Vector2d(1, 0));
    auto& n2 = geoMesh.CoordinateNodes.Add(Eigen::Vector2d(0, 1));
    auto& n3 = geoMesh.CoordinateNodes.Add(Eigen::Vector2d(1, 1));

    auto& interpolation = mesh.CreateInterpolation(NuTo::InterpolationTriangleLinear());
    auto& cElm0 = geoMesh.Elements.Add({{n0, n1, n2}, interpolation});
    auto& cElm1 = geoMesh.Elements.Add({{n1, n3, n2}, interpolation});

    mesh.Elements.Add(cElm0);
    mesh.Elements.Add(cElm1);

    int expectedNumCoordinateNodes = 4;
    BOOST_CHECK_EQUAL(geoMesh.CoordinateNodes.Size(), expectedNumCoordinateNodes);


    // add linear dof type

    NuTo::DofType dof0("linear", 1);
    NuTo::AddDofInterpolation(&mesh, dof0, interpolation);

    int expectedNumDof0Nodes = expectedNumCoordinateNodes; // same interpolation

    BOOST_CHECK_EQUAL(mesh.Nodes.Size(), expectedNumDof0Nodes);
    BOOST_CHECK_EQUAL(geoMesh.CoordinateNodes.Size(), expectedNumCoordinateNodes);
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
    BOOST_CHECK_EQUAL(mesh.Nodes.Size(), expectedNumDof0Nodes + expectedNumDof1Nodes);
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
    NuTo::GeometryMeshFem geoMesh;
    NuTo::MeshFem mesh(geoMesh);
    auto& n0 = geoMesh.CoordinateNodes.Add(Eigen::Vector2d(0, 0));
    auto& n1 = geoMesh.CoordinateNodes.Add(Eigen::Vector2d(1, 0));
    auto& n2 = geoMesh.CoordinateNodes.Add(Eigen::Vector2d(0, 1));
    auto& n3 = geoMesh.CoordinateNodes.Add(Eigen::Vector2d(1, 1));
    auto& n4 = geoMesh.CoordinateNodes.Add(Eigen::Vector2d(2, 0));
    auto& n5 = geoMesh.CoordinateNodes.Add(Eigen::Vector2d(2, 1));

    auto& interpolationTriangle = mesh.CreateInterpolation(NuTo::InterpolationTriangleLinear());
    auto& interpolationQuad = mesh.CreateInterpolation(NuTo::InterpolationQuadLinear());
    auto& cElm0 = geoMesh.Elements.Add({{n0, n1, n2}, interpolationTriangle});
    auto& cElm1 = geoMesh.Elements.Add({{n1, n4, n5, n3}, interpolationQuad});
    mesh.Elements.Add(cElm0);
    mesh.Elements.Add(cElm1);

    int expectedNumCoordinateNodes = 6;
    BOOST_CHECK_EQUAL(geoMesh.CoordinateNodes.Size(), expectedNumCoordinateNodes);


    // add linear dof type

    NuTo::DofType dof0("isoparametric", 1);
    NuTo::AddDofInterpolation(&mesh, dof0);

    int expectedNumDof0Nodes = expectedNumCoordinateNodes; // same interpolation

    BOOST_CHECK_EQUAL(mesh.Nodes.Size(), expectedNumDof0Nodes);
    BOOST_CHECK_NO_THROW(mesh.NodeAtCoordinate(Eigen::Vector2d(0, 0), dof0));
}

BOOST_AUTO_TEST_CASE(MeshNodesTotalDof)
{
    NuTo::GeometryMeshFem geoMesh;
    NuTo::MeshFem mesh(geoMesh);
    auto& n0 = geoMesh.CoordinateNodes.Add(Eigen::Vector2d(0, 0));
    auto& n1 = geoMesh.CoordinateNodes.Add(Eigen::Vector2d(1, 0));
    auto& n2 = geoMesh.CoordinateNodes.Add(Eigen::Vector2d(0, 1));

    auto& interpolationTriangle = mesh.CreateInterpolation(NuTo::InterpolationTriangleLinear());
    auto& interpolationTruss = mesh.CreateInterpolation(NuTo::InterpolationTrussLinear());

    // Add coordinate elements
    auto& cElm0 = geoMesh.Elements.Add({{n0, n1}, interpolationTruss});
    auto& line1 = mesh.Elements.Add(cElm0);

    auto& cElm1 = geoMesh.Elements.Add({{n0, n1, n2}, interpolationTriangle});
    auto& cElm2 = geoMesh.Elements.Add({{n1, n2}, interpolationTruss});
    auto& cElm3 = geoMesh.Elements.Add({{n2, n0}, interpolationTruss});

    mesh.Elements.Add(cElm1);
    mesh.Elements.Add(cElm2);
    mesh.Elements.Add(cElm3);

    // Add a dof element
    NuTo::DofType dof1("dof1", 1);

    auto& nD0 = mesh.Nodes.Add(1.2);
    auto& nD1 = mesh.Nodes.Add(2.2);

    line1.AddDofElement(dof1, {{nD0, nD1}, interpolationTruss});

    BOOST_CHECK_EQUAL(mesh.NodesTotal(dof1).Size(), 2);
}

BOOST_AUTO_TEST_CASE(PartialAddDofConvert)
{
    NuTo::GeometryMeshFem geoMesh;
    NuTo::MeshFem mesh(geoMesh);
    auto& n0 = geoMesh.CoordinateNodes.Add(Eigen::Vector2d(0, 0));
    auto& n1 = geoMesh.CoordinateNodes.Add(Eigen::Vector2d(1, 0));
    auto& n2 = geoMesh.CoordinateNodes.Add(Eigen::Vector2d(0, 1));

    auto& interpolationTriangle = mesh.CreateInterpolation(NuTo::InterpolationTriangleLinear());
    auto& interpolationTruss = mesh.CreateInterpolation(NuTo::InterpolationTrussLinear());

    // Add coordinate elements
    auto& cElm1 = geoMesh.Elements.Add({{n0, n1, n2}, interpolationTriangle});
    auto& cElm2 = geoMesh.Elements.Add({{n0, n1}, interpolationTruss});
    auto& cElm3 = geoMesh.Elements.Add({{n1, n2}, interpolationTruss});
    auto& cElm4 = geoMesh.Elements.Add({{n2, n0}, interpolationTruss});

    auto& tri = mesh.Elements.Add(cElm1);
    auto& line1 = mesh.Elements.Add(cElm2);
    auto& line2 = mesh.Elements.Add(cElm3);
    auto& line3 = mesh.Elements.Add(cElm4);

    // Add dof element in steps
    NuTo::DofType dof1("dof1", 1);

    NuTo::AddDofInterpolation(&mesh, dof1, {line1}, interpolationTruss);
    BOOST_CHECK_EQUAL(mesh.NodesTotal(dof1).Size(), 2);

    NuTo::AddDofInterpolation(&mesh, dof1, {line2, line3}, interpolationTruss);
    BOOST_CHECK_EQUAL(mesh.NodesTotal(dof1).Size(), 3);

    NuTo::AddDofInterpolation(&mesh, dof1, {tri}, interpolationTriangle);
    BOOST_CHECK_EQUAL(mesh.NodesTotal(dof1).Size(), 3);

    // throw if the interpolation doesn't match the shape of the underlying elements
    BOOST_CHECK_THROW(NuTo::AddDofInterpolation(&mesh, dof1, {tri}, interpolationTruss), NuTo::Exception);
}
