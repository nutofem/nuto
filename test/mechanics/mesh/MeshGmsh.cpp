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
    BOOST_CHECK_NO_THROW("quad1.msh");
}

BOOST_AUTO_TEST_CASE(QuadSerendipity)
{
    BOOST_CHECK_NO_THROW("quadSerendipity2.msh");
}

BOOST_AUTO_TEST_CASE(TriangleLinear)
{
    BOOST_CHECK_NO_THROW("triangle1.msh");
}

BOOST_AUTO_TEST_CASE(TriangleQuadratic)
{
    BOOST_CHECK_NO_THROW("triangle2.msh");
}

BOOST_AUTO_TEST_CASE(LineLinear)
{
    BOOST_CHECK_NO_THROW("line1.msh");
}

BOOST_AUTO_TEST_CASE(TetrahedronLinear)
{
    BOOST_CHECK_NO_THROW("tet1.msh");
}

BOOST_AUTO_TEST_CASE(BrickLinear)
{
    BOOST_CHECK_NO_THROW("hex1.msh");
}

BOOST_AUTO_TEST_CASE(PrismLinear)
{
    BOOST_CHECK_NO_THROW("prism1.msh");
}

BOOST_AUTO_TEST_CASE(PyramidLinear)
{
    BOOST_CHECK_NO_THROW("pyramid1.msh");
}
