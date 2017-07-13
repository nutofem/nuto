#include "DamageLawHelper.h"
#include "mechanics/constitutive/damageLaws/DamageLawLinear.h"

BOOST_AUTO_TEST_CASE(LinearDerivative)
{
    DamageLawHelper::CheckDerivatives(*NuTo::Constitutive::DamageLawLinear::Create(1e-4, 0.4, 0.5));
    DamageLawHelper::CheckDerivatives(*NuTo::Constitutive::DamageLawLinear::Create(1e-4, 0.4, 1.));
}

BOOST_AUTO_TEST_CASE(LinearZerosKappaC)
{
    auto law = NuTo::Constitutive::DamageLawLinear::Create(1e-4, 0.4, 1);
    for (double kappa : {0.4, 0.5, 0.6})
    {
        BOOST_CHECK_CLOSE(law->CalculateDamage(kappa), 1., 1.e-10);
        BOOST_CHECK_SMALL(law->CalculateDerivative(kappa), 1.e-10);
    }
}

BOOST_AUTO_TEST_CASE(LinearDamageMax)
{
    for (double damageMax : {0.5, 0.9, 1.})
    {
        auto law = NuTo::Constitutive::DamageLawLinear::Create(1e-4, 0.4, damageMax);
        for (double kappa = 0; kappa < 0.5; kappa += 1.e-4)
            BOOST_CHECK_LE(law->CalculateDamage(kappa), damageMax);
    }
}
