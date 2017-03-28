#pragma once

#include <memory>
#include "mechanics/constitutive/damageLaws/DamageLaw.h"

namespace NuTo
{
namespace Constitutive
{

class DamageLawLinear : public DamageLaw
{
public:
    static std::shared_ptr<DamageLaw> Create(double kappa0, double kappaC, double damageMax)
    {
        return std::shared_ptr<DamageLaw>(new DamageLawLinear(kappa0, kappaC, damageMax));
    }

protected:
    DamageLawLinear(double kappa0, double kappaC, double damageMax)
        : DamageLaw(kappa0)
        , mKappaC(kappaC)
        , mDamageMax(damageMax)
    {
    }

    //! @brief protected virtual method for the damage calculation
    //!        the case kappa <= kappa0 is already covered in the public interface
    //! @param kappa history variable
    //! @return damage
    double Damage(const double kappa) const
    {
        double damage = mKappaC / kappa * (kappa - mKappa0) / (mKappaC - mKappa0);
        return std::min(damage, mDamageMax);
    }

    //! @brief protected virtual method for the damage derivative calculation
    //!        the case kappa <= kappa0 is already covered in the public interface
    //! @param kappa history variable
    //! @return damage derivative
    double Derivative(const double kappa) const
    {
        if (Damage(kappa) < mDamageMax)
            return mKappaC * mKappa0 / (kappa * kappa * (mKappaC - mKappa0));
        return 0;
    }

private:
    const double mKappaC;
    const double mDamageMax;
};

} /* Constitutive */
} /* NuTo */
