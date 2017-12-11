#include "BoostUnitTest.h"
#include "TypeTraits.h"
#include "mechanics/constitutive/EngineeringStrain.h"

BOOST_AUTO_TEST_CASE(EngineeringStrainCopyMove)
{
    NuTo::Test::Copy<NuTo::EngineeringStrain<3>>();
    NuTo::Test::Move<NuTo::EngineeringStrain<3>>();
}
