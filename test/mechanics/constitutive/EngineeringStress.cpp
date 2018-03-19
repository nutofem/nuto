#include "BoostUnitTest.h"
#include "TypeTraits.h"
#include "nuto/mechanics/constitutive/EngineeringStress.h"

BOOST_AUTO_TEST_CASE(EngineeringStrainCopyMove)
{
    NuTo::Test::Copy<NuTo::EngineeringStress<3>>();
    NuTo::Test::Move<NuTo::EngineeringStress<3>>();
}
