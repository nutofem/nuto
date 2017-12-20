#include "mechanics/mesh/MeshGmsh.h"
#include "BoostUnitTest.h"

using namespace NuTo;

BOOST_AUTO_TEST_CASE(PhysicalGroups)
{
    MeshGmsh gmsh("quad.msh");
    BOOST_CHECK_NO_THROW(gmsh.GetPhysicalGroup("domain"));
    BOOST_CHECK_THROW(gmsh.GetPhysicalGroup("are_you_there?"), Exception);
}
