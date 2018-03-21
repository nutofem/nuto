#include "nuto/mechanics/mesh/UnitMeshFem.h"
#include "nuto/mechanics/cell/Jacobian.h"
#include "BoostUnitTest.h"

void CheckJacobians(NuTo::MeshFem& mesh)
{
    int dim = mesh.Elements[0].CoordinateElement().GetNode(0).GetValues().rows();
    Eigen::VectorXd ip = Eigen::VectorXd::Zero(dim);
    for (auto& element : mesh.Elements)
    {
        auto d_dxi = element.CoordinateElement().GetDerivativeShapeFunctions(ip);
        auto x = element.CoordinateElement().ExtractNodeValues();
        auto J = NuTo::Jacobian(x, d_dxi, dim);
        BOOST_CHECK_GT(J.Det(), 0.);
    }
}

void Check2DMesh(NuTo::MeshFem& mesh)
{
    BOOST_CHECK_EQUAL(mesh.Nodes.Size(), 3 * 8);

    BOOST_CHECK_NO_THROW(mesh.NodeAtCoordinate(Eigen::Vector2d(0., 0.)));
    BOOST_CHECK_NO_THROW(mesh.NodeAtCoordinate(Eigen::Vector2d(1., 1.)));

    BOOST_CHECK_NO_THROW(mesh.NodeAtCoordinate(Eigen::Vector2d(1. / 2., 1. / 7.)));
    BOOST_CHECK_NO_THROW(mesh.NodeAtCoordinate(Eigen::Vector2d(1. / 2., 5. / 7.)));

    CheckJacobians(mesh);

    auto f = [](Eigen::VectorXd coords) // transforms mesh to (42,4) -- (44, 11)
    {
        Eigen::VectorXd newCoords(2);
        newCoords[0] = 42 + coords[0] * 2;
        newCoords[1] = 4 + coords[1] * 7;
        return newCoords;
    };

    NuTo::MeshFem transformedMesh = NuTo::UnitMeshFem::Transform(std::move(mesh), f);

    BOOST_CHECK_NO_THROW(transformedMesh.NodeAtCoordinate(Eigen::Vector2d(42., 4.)));
    BOOST_CHECK_NO_THROW(transformedMesh.NodeAtCoordinate(Eigen::Vector2d(44., 11.)));

    BOOST_CHECK_NO_THROW(transformedMesh.NodeAtCoordinate(Eigen::Vector2d(43., 5.)));
    BOOST_CHECK_NO_THROW(transformedMesh.NodeAtCoordinate(Eigen::Vector2d(43., 9.)));
}

BOOST_AUTO_TEST_CASE(MeshTrusses)
{
    constexpr int numElements = 15;
    auto mesh = NuTo::UnitMeshFem::CreateLines(numElements);
    BOOST_CHECK_EQUAL(mesh.Elements.Size(), numElements);
    BOOST_CHECK_EQUAL(mesh.Nodes.Size(), numElements + 1);

    auto IsWholeNumber = [](double d, double eps = 1.e-12) { return std::abs(d - std::floor(d)) < eps; };

    for (const auto& node : mesh.Nodes)
    {
        BOOST_CHECK(IsWholeNumber(node.GetValues()[0] * numElements));
        BOOST_CHECK_LE(node.GetValues()[0], 1.0);
        BOOST_CHECK_GE(node.GetValues()[0], 0.0);
    }

    for (const auto& element : mesh.Elements)
    {
        BOOST_CHECK_LT(element.CoordinateElement().GetNode(0).GetValues()[0],
                       element.CoordinateElement().GetNode(1).GetValues()[0]);
    }
}

BOOST_AUTO_TEST_CASE(MeshQuad)
{
    auto mesh = NuTo::UnitMeshFem::CreateQuads(2, 7);
    BOOST_CHECK_EQUAL(mesh.Elements.Size(), 2 * 7);
    Check2DMesh(mesh);
}

BOOST_AUTO_TEST_CASE(MeshTriangle)
{
    auto mesh = NuTo::UnitMeshFem::CreateTriangles(2, 7);
    BOOST_CHECK_EQUAL(mesh.Elements.Size(), 2 * 7 * 2);
    Check2DMesh(mesh);
}

BOOST_AUTO_TEST_CASE(MeshBrick)
{
    auto mesh = NuTo::UnitMeshFem::CreateBricks(2, 7, 3);
    BOOST_CHECK_EQUAL(mesh.Elements.Size(), 2 * 7 * 3);
    BOOST_CHECK_EQUAL(mesh.Nodes.Size(), 3 * 8 * 4);
    BOOST_CHECK_NO_THROW(mesh.NodeAtCoordinate(Eigen::Vector3d(0, 0, 0)));
    BOOST_CHECK_NO_THROW(mesh.NodeAtCoordinate(Eigen::Vector3d(1, 1, 1)));
    CheckJacobians(mesh);
}

BOOST_AUTO_TEST_CASE(MeshValidAfterTransform)
{
    auto mesh = NuTo::UnitMeshFem::CreateQuads(1, 1);
    Eigen::VectorXd expected(8);
    expected << 0, 0, 1, 0, 1, 1, 0, 1;

    auto& coordinateElement = mesh.Elements[0].CoordinateElement();
    BoostUnitTest::CheckEigenMatrix(coordinateElement.ExtractNodeValues(), expected);

    auto f = [](Eigen::VectorXd coords) { return Eigen::Vector2d(coords[0] * 4, coords[1] * 42); };

    NuTo::MeshFem transformedMesh = NuTo::UnitMeshFem::Transform(std::move(mesh), f);
    auto& transformedCoordinateElement = transformedMesh.Elements[0].CoordinateElement();
    expected << 0, 0, 4, 0, 4, 42, 0, 42;
    BoostUnitTest::CheckEigenMatrix(transformedCoordinateElement.ExtractNodeValues(), expected);

    transformedMesh.Nodes[0].SetValue(0, 6174);
    expected << 6174, 0, 4, 0, 4, 42, 0, 42;
    BoostUnitTest::CheckEigenMatrix(transformedCoordinateElement.ExtractNodeValues(), expected);
}

// This test is related to our github issue #148. Visit github to read about the details
BOOST_AUTO_TEST_CASE(MeshMovabilityError)
{
    NuTo::MeshFem mesh = NuTo::UnitMeshFem::CreateLines(1);
    {
        NuTo::MeshFem tempMesh = NuTo::UnitMeshFem::CreateLines(1);
        mesh = std::move(tempMesh);
    }
    auto& coordinateElement = mesh.Elements[0].CoordinateElement();
    coordinateElement.GetDofDimension();
}
