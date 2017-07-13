#include "DamageLawHelper.h"
#include "mechanics/constitutive/damageLaws/DamageLawHermite.h"

BOOST_AUTO_TEST_CASE(HermiteDerivative)
{
    DamageLawHelper::CheckDerivatives(*NuTo::Constitutive::DamageLawHermite::Create(1e-4, 0.4));
}

BOOST_AUTO_TEST_CASE(HermiteZerosKappaC)
{
    auto law = NuTo::Constitutive::DamageLawHermite::Create(1e-4, 0.4);
    for (double kappa : {0.4, 0.5, 0.6})
    {
        BOOST_CHECK_CLOSE(law->CalculateDamage(kappa), 1., 1.e-10);
        BOOST_CHECK_SMALL(law->CalculateDerivative(kappa), 1.e-10);
    }
}
