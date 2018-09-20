#include "nuto/mechanics/mesh/MeshGmsh.h"
#include "BoostUnitTest.h"
#include "nuto/mechanics/cell/Jacobian.h"

using namespace NuTo;

BOOST_AUTO_TEST_CASE(PhysicalGroups)
{
    MeshGmsh gmsh("meshes/quad.msh");

    BOOST_CHECK_NO_THROW(gmsh.GetPhysicalGroup(1));
    BOOST_CHECK_THROW(gmsh.GetPhysicalGroup(-42), Exception);

    BOOST_CHECK_NO_THROW(gmsh.GetPhysicalGroup("domain"));
    BOOST_CHECK_THROW(gmsh.GetPhysicalGroup("are_you_there?"), Exception);
}

BOOST_AUTO_TEST_CASE(NonContiguousNodeNumbering)
{
    BOOST_CHECK_NO_THROW(MeshGmsh("meshes/quadNoncontiguous.msh"));
}

void CheckMesh(std::string meshFile, int numNodesExpected)
{
    BOOST_CHECK_NO_THROW(MeshGmsh{meshFile});
    MeshGmsh m(meshFile);
    auto& meshFem = m.GetMeshFEM();
    BOOST_CHECK_EQUAL(meshFem.NumNodes(), numNodesExpected);
    const CoordinateElementFem& element = *meshFem.GetElements().begin();
    const auto& interpolation = element.Interpolation();
    BOOST_CHECK_EQUAL(interpolation.GetNumNodes(), numNodesExpected);

    // check for consistent (same sign for all IP) Det(Jacobian)
    bool isPositive = true;
    for (int iNode = 0; iNode < interpolation.GetNumNodes(); ++iNode)
    {
        auto coord = interpolation.GetLocalCoords(iNode);
        Jacobian j(element.ExtractCoordinates(), interpolation.GetDerivativeShapeFunctions(coord));
        bool isPositiveForINode = j.Det() > -1.e-10;
        if (iNode == 0)
            isPositive = isPositiveForINode;
        else
            BOOST_CHECK(isPositive == isPositiveForINode);
    }
}

BOOST_AUTO_TEST_CASE(BinaryImport)
{
    MeshGmsh m("meshes/binary.msh");
    auto& meshFem = m.GetMeshFEM();
    BOOST_CHECK_EQUAL(meshFem.NumNodes(), 8);
    BOOST_CHECK_EQUAL(meshFem.GetElements().Size(), 7);
    BOOST_CHECK_EQUAL(meshFem.GetElements().begin()->Interpolation().GetNumNodes(), 3);
}

BOOST_AUTO_TEST_CASE(QuadLinear)
{
    CheckMesh("meshes/quad1.msh", 4);
}


BOOST_AUTO_TEST_CASE(QuadSerendipity)
{
    CheckMesh("meshes/quadSerendipity2.msh", 8);
}

BOOST_AUTO_TEST_CASE(TriangleLinear)
{
    CheckMesh("meshes/triangle1.msh", 3);
}

BOOST_AUTO_TEST_CASE(TriangleQuadratic)
{
    CheckMesh("meshes/triangle2.msh", 6);
}

BOOST_AUTO_TEST_CASE(LineLinear)
{
    CheckMesh("meshes/line1.msh", 2);
}

BOOST_AUTO_TEST_CASE(LineQuadratic)
{
    CheckMesh("meshes/line2.msh", 3);
}

BOOST_AUTO_TEST_CASE(TetrahedronLinear)
{
    CheckMesh("meshes/tet1.msh", 4);
}

BOOST_AUTO_TEST_CASE(BrickLinear)
{
    CheckMesh("meshes/hex1.msh", 8);
}

BOOST_AUTO_TEST_CASE(PrismLinear)
{
    CheckMesh("meshes/prism1.msh", 6);
}

BOOST_AUTO_TEST_CASE(PyramidLinear)
{
    CheckMesh("meshes/pyramid1.msh", 5);
}

BOOST_AUTO_TEST_CASE(QuadLobatto2)
{
    CheckMesh("meshes/quadL2.msh", 9);
}

BOOST_AUTO_TEST_CASE(BrickLobatto2)
{
    CheckMesh("meshes/brickL2.msh", 27);
}

BOOST_AUTO_TEST_CASE(InvertedTriangle)
{
    // triangle mesh with negative Jacobian determinant
    BOOST_CHECK_THROW(MeshGmsh("meshes/invertedTriangle.msh"), Exception);
}
