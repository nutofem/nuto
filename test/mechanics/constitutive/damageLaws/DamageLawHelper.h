#pragma once
#include "BoostUnitTest.h"
#include "mechanics/constitutive/damageLaws/DamageLaw.h"

namespace DamageLawHelper
{
//! @brief checks the derivatives of the damage law via central differences
//! @param law damage law
void CheckDerivatives(const NuTo::Constitutive::DamageLaw& law, double kappaEnd = .5)
{
    const double delta = 1.e-8;
    const double kappa0 = law.GetKappa0();
    for (double kappa = delta; kappa < kappaEnd; kappa += kappa0/5.)
    {
        if (kappa < kappa0)
        {
            BOOST_CHECK_SMALL(law.CalculateDamage(kappa), 1.e-10);
            BOOST_CHECK_SMALL(law.CalculateDerivative(kappa), 1.e-10);
        }

        const double derivative = law.CalculateDerivative(kappa);
        const double derivativeCDF =
                (law.CalculateDamage(kappa + .5 * delta) - law.CalculateDamage(kappa - .5 * delta)) / delta;
        BOOST_CHECK_SMALL(derivative - derivativeCDF, 1.e-4);
    }
}
} /* DamageLawHelper */
