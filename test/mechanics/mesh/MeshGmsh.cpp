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

BOOST_AUTO_TEST_CASE(QuadLinear)
{
    MeshGmsh gmsh("quad1.msh");
    auto g = gmsh.GetPhysicalGroup(0);
}

BOOST_AUTO_TEST_CASE(QuadSerendipity)
{
    MeshGmsh gmsh("quadSerendipity2.msh");
    auto g = gmsh.GetPhysicalGroup(0);
}

BOOST_AUTO_TEST_CASE(TriangleLinear)
{
    MeshGmsh gmsh("triangle1.msh");
    auto g = gmsh.GetPhysicalGroup(0);
}

BOOST_AUTO_TEST_CASE(TriangleQuadratic)
{
    MeshGmsh gmsh("triangle2.msh");
    auto g = gmsh.GetPhysicalGroup(0);
}

BOOST_AUTO_TEST_CASE(LineLinear)
{
    MeshGmsh gmsh("line1.msh");
    auto g = gmsh.GetPhysicalGroup(0);
}

BOOST_AUTO_TEST_CASE(TetrahedronLinear)
{
    MeshGmsh gmsh("tet1.msh");
    auto g = gmsh.GetPhysicalGroup(0);
}

BOOST_AUTO_TEST_CASE(BrickLinear)
{
    MeshGmsh gmsh("hex1.msh");
    auto g = gmsh.GetPhysicalGroup(0);
}

BOOST_AUTO_TEST_CASE(PrismLinear)
{
    MeshGmsh gmsh("prism1.msh");
    auto g = gmsh.GetPhysicalGroup(0);
}

BOOST_AUTO_TEST_CASE(PyramidLinear)
{
    MeshGmsh gmsh("pyramid1.msh");
    auto g = gmsh.GetPhysicalGroup(0);
}
