#pragma once
#include "BoostUnitTest.h"
#include "mechanics/constitutive/damageLaws/DamageLaw.h"

namespace DamageLawHelper
{
//! @brief checks the derivatives of the damage law via central differences
//! @param law damage law
//! @param kappa0 initial damage threshold
//! @param kappaEnd kappa range
void CheckDerivatives(const NuTo::Constitutive::DamageLaw& law, double kappa0, double kappaEnd = .5)
{
    const double delta = 1.e-8;
    for (double kappa = delta; kappa < kappaEnd; kappa += kappa0 / 5.)
    {
        if (kappa < kappa0)
        {
            BOOST_CHECK_SMALL(law.Damage(kappa), 1.e-10);
            BOOST_CHECK_SMALL(law.Derivative(kappa), 1.e-10);
        }
        const double derivative = law.Derivative(kappa);
        const double derivativeCDF = (law.Damage(kappa + .5 * delta) - law.Damage(kappa - .5 * delta)) / delta;
        BOOST_CHECK_SMALL(derivative - derivativeCDF, 1.e-4);
    }
}
} /* DamageLawHelper */
