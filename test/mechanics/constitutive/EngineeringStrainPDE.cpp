#include "BoostUnitTest.h"
#include "TypeTraits.h"
#include "mechanics/constitutive/EngineeringStrainPDE.h"

BOOST_AUTO_TEST_CASE(EngineeringStrainCopyMove)
{
    NuTo::Test::Copy<NuTo::EngineeringStrainPDE<3>>();
    NuTo::Test::Move<NuTo::EngineeringStrainPDE<3>>();
}
