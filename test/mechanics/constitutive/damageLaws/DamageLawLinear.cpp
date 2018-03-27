#include "DamageLawHelper.h"
#include "nuto/mechanics/constitutive/damageLaws/DamageLawLinear.h"

BOOST_AUTO_TEST_CASE(LinearDerivative)
{
    DamageLawHelper::CheckDerivatives(NuTo::Constitutive::DamageLawLinear(1e-4, 0.4, 0.5), 1e-4);
    DamageLawHelper::CheckDerivatives(NuTo::Constitutive::DamageLawLinear(1e-4, 0.4, 1.), 1e-4);
}

BOOST_AUTO_TEST_CASE(LinearZerosKappaC)
{
    NuTo::Constitutive::DamageLawLinear law(1e-4, 0.4, 1);
    for (double kappa : {0.4, 0.5, 0.6})
    {
        BOOST_CHECK_CLOSE(law.Damage(kappa), 1., 1.e-10);
        BOOST_CHECK_SMALL(law.Derivative(kappa), 1.e-10);
    }
}

BOOST_AUTO_TEST_CASE(LinearDamageMax)
{
    for (double damageMax : {0.5, 0.9, 1.})
    {
        NuTo::Constitutive::DamageLawLinear law(1e-4, 0.4, damageMax);
        for (double kappa = 0; kappa < 0.5; kappa += 1.e-4)
            BOOST_CHECK_LE(law.Damage(kappa), damageMax);
    }
}

BOOST_AUTO_TEST_CASE(LinearDamageMaterial)
{
    NuTo::Material::Softening m = NuTo::Material::DefaultConcrete();
    m.fMin = 0;
    NuTo::Constitutive::DamageLawLinear law(m);
    DamageLawHelper::CheckFractureEnergy(law, m);
}
