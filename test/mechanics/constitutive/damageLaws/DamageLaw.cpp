#include "BoostUnitTest.h"
#include "mechanics/constitutive/damageLaws/DamageLaw.h"

class MockLaw : public NuTo::Constitutive::DamageLaw
{
public:
    MockLaw()
        : DamageLaw(1.e-4)
    {
    }

protected:
    double Damage(const double) const override
    {
        return 42;
    }
    double Derivative(const double) const override
    {
        return 42;
    }
};

//! @brief checks if the damage and the derivative is actually zero below kappa0
//! @param law damage law
BOOST_AUTO_TEST_CASE(Zeros)
{
    MockLaw law;
    const double kappa0 = law.GetKappa0();
    BOOST_CHECK_CLOSE(kappa0, 1.e-4, 1.e-10);

    for (double kappa : {0., 0.5 * kappa0, kappa0})
    {
        BOOST_CHECK_SMALL(law.CalculateDamage(kappa), 1.e-10);
        BOOST_CHECK_SMALL(law.CalculateDerivative(kappa), 1.e-10);
    }
}
