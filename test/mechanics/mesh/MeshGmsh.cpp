#include "mechanics/mesh/MeshGmsh.h"
#include "BoostUnitTest.h"

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
    BOOST_CHECK_EQUAL(meshFem.Nodes.Size(), numNodesExpected);
    BOOST_CHECK_EQUAL(meshFem.Elements.begin()->CoordinateElement().Interpolation().GetNumNodes(), numNodesExpected);
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
