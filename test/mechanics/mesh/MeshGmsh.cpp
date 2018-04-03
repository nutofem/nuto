#include "nuto/mechanics/mesh/MeshGmsh.h"
#include "BoostUnitTest.h"
#include "nuto/mechanics/cell/Jacobian.h"

using namespace NuTo;

BOOST_AUTO_TEST_CASE(PhysicalGroups)
{
    MeshGmsh gmsh("quad.msh");

    BOOST_CHECK_NO_THROW(gmsh.GetPhysicalGroup(1));
    BOOST_CHECK_THROW(gmsh.GetPhysicalGroup(-42), Exception);

    BOOST_CHECK_NO_THROW(gmsh.GetPhysicalGroup("domain"));
    BOOST_CHECK_THROW(gmsh.GetPhysicalGroup("are_you_there?"), Exception);
}

BOOST_AUTO_TEST_CASE(NonContiguousNodeNumbering)
{
    BOOST_CHECK_NO_THROW(MeshGmsh("quadNoncontiguous.msh"));
}

void CheckMesh(std::string meshFile, int numNodesExpected)
{
    BOOST_CHECK_NO_THROW(MeshGmsh{meshFile});
    MeshGmsh m(meshFile);
    auto& meshFem = m.GetMeshFEM();
    BOOST_CHECK_EQUAL(meshFem.CoordinateNodes.Size(), numNodesExpected);
    const CoordinateElementFem& element = meshFem.Elements.begin()->CoordinateElement();
    const auto& interpolation = element.Interpolation();
    BOOST_CHECK_EQUAL(interpolation.GetNumNodes(), numNodesExpected);

    // check for consistent (same sign for all IP) Det(Jacobian)
    bool isPositive = true;
    for (int iNode = 0; iNode < interpolation.GetNumNodes(); ++iNode)
    {
        auto coord = interpolation.GetLocalCoords(iNode);
        int globalDimension = coord.rows();
        Jacobian j(element.ExtractNodeValues(), interpolation.GetDerivativeShapeFunctions(coord), globalDimension);
        bool isPositiveForINode = j.Det() > -1.e-10;
        if (iNode == 0)
            isPositive = isPositiveForINode;
        else
            BOOST_CHECK(isPositive == isPositiveForINode);
    }
}

BOOST_AUTO_TEST_CASE(BinaryImport)
{
    MeshGmsh m("binary.msh");
    auto& meshFem = m.GetMeshFEM();
    BOOST_CHECK_EQUAL(meshFem.CoordinateNodes.Size(), 8);
    BOOST_CHECK_EQUAL(meshFem.Elements.Size(), 7);
    BOOST_CHECK_EQUAL(meshFem.Elements.begin()->CoordinateElement().Interpolation().GetNumNodes(), 3);
}

BOOST_AUTO_TEST_CASE(QuadLinear)
{
    CheckMesh("quad1.msh", 4);
}


BOOST_AUTO_TEST_CASE(QuadSerendipity)
{
    CheckMesh("quadSerendipity2.msh", 8);
}

BOOST_AUTO_TEST_CASE(TriangleLinear)
{
    CheckMesh("triangle1.msh", 3);
}

BOOST_AUTO_TEST_CASE(TriangleQuadratic)
{
    CheckMesh("triangle2.msh", 6);
}

BOOST_AUTO_TEST_CASE(LineLinear)
{
    CheckMesh("line1.msh", 2);
}

BOOST_AUTO_TEST_CASE(LineQuadratic)
{
    CheckMesh("line2.msh", 3);
}

BOOST_AUTO_TEST_CASE(TetrahedronLinear)
{
    CheckMesh("tet1.msh", 4);
}

BOOST_AUTO_TEST_CASE(BrickLinear)
{
    CheckMesh("hex1.msh", 8);
}

BOOST_AUTO_TEST_CASE(PrismLinear)
{
    CheckMesh("prism1.msh", 6);
}

BOOST_AUTO_TEST_CASE(PyramidLinear)
{
    CheckMesh("pyramid1.msh", 5);
}

BOOST_AUTO_TEST_CASE(QuadLobatto2)
{
    CheckMesh("quadL2.msh", 9);
}

BOOST_AUTO_TEST_CASE(BrickLobatto2)
{
    CheckMesh("brickL2.msh", 27);
}
