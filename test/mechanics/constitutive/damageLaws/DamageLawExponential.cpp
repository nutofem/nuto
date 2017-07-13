#include "DamageLawHelper.h"
#include "mechanics/constitutive/damageLaws/DamageLawExponential.h"

BOOST_AUTO_TEST_CASE(ExponentialDerivative)
{
    DamageLawHelper::CheckDerivatives(*NuTo::Constitutive::DamageLawExponential::Create(1e-4, 350, 0.9));
}

BOOST_AUTO_TEST_CASE(ExponentialResLoad)
{
    const double kappa0 = 1.e-4;
    const double alpha = 0.9;
    const double E = 1.e6;
    const double sigma_infty = E * kappa0 * (1. - alpha);
    BOOST_TEST_MESSAGE("sigma" << sigma_infty);
    auto law = NuTo::Constitutive::DamageLawExponential::Create(1.e-4, 350, alpha);
    for (double kappa : {1., 2., 10., 100., 1000.})
    {
        const double pseudoStress = (1. - law->CalculateDamage(kappa)) * E * kappa;
        BOOST_CHECK_CLOSE(pseudoStress, sigma_infty, 1.e-3);
    }
}
