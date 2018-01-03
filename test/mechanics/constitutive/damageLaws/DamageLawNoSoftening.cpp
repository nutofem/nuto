#include "DamageLawHelper.h"
#include "mechanics/constitutive/damageLaws/DamageLawNoSoftening.h"

BOOST_AUTO_TEST_CASE(LinearDerivative)
{
    DamageLawHelper::CheckDerivatives(NuTo::Constitutive::DamageLawNoSoftening(1e-4), 1.e-4);
}
