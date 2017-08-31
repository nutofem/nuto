#include "BoostUnitTest.h"
#include "TypeTraits.h"
#include "mechanics/constitutive/EngineeringStressPDE.h"

BOOST_AUTO_TEST_CASE(EngineeringStrainCopyMove)
{
    NuTo::Test::Copy<NuTo::EngineeringStressPDE<3>>();
    NuTo::Test::Move<NuTo::EngineeringStressPDE<3>>();
}
