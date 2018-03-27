#pragma once
#include "BoostUnitTest.h"
#include "nuto/mechanics/constitutive/damageLaws/DamageLaw.h"
#include "nuto/mechanics/constitutive/damageLaws/SofteningMaterial.h"

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

void CheckFractureEnergy(const NuTo::Constitutive::DamageLaw& law, NuTo::Material::Softening m)
{
    double k0 = m.ft / m.E;
    double dk = k0 / 10.;
    double gf = 0;
    for (double k = k0; k < 10000 * k0; k += dk)
    {
        double sigma_a = (1. - law.Damage(k)) * m.E * k;
        double sigma_b = (1. - law.Damage(k + dk)) * m.E * (k + dk);
        gf += 0.5 * dk * (sigma_a + sigma_b);
    }
    BOOST_CHECK_CLOSE(gf, m.gf, 1.e-5);
}
} /* DamageLawHelper */
