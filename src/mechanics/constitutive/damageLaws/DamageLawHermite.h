#pragma once

#include <memory>
#include <cmath>
#include "mechanics/constitutive/damageLaws/DamageLaw.h"

namespace NuTo
{
namespace Constitutive
{

class DamageLawHermite : public DamageLaw
{
public:
    //! @brief Create method
    //! @param kappa0 initial kappa
    //! @param alpha max damage
    //! @param beta represents tensile strength / fracture energy
    static std::shared_ptr<DamageLaw> Create(double kappa0, double kappaC)
    {
        return std::shared_ptr<DamageLaw>(new DamageLawHermite(kappa0, kappaC));
    }

protected:
    DamageLawHermite(double kappa0, double kappaC)
        : DamageLaw(kappa0)
        , mKappaC(kappaC)
    {
    }

    //! @brief protected virtual method for the damage calculation
    //!        the case kappa <= kappa0 is already covered in the public interface
    //! @param kappa history variable
    //! @return damage
    double Damage(const double kappa) const
    {
        if (kappa >= mKappaC)
            return 1;
        const double kappaScaled = GetKappaScaled(kappa);
        return 1 - mKappa0 / kappa * (2 * kappaScaled * kappaScaled * kappaScaled - 3 * kappaScaled * kappaScaled + 1);
    }

    //! @brief protected virtual method for the damage derivative calculation
    //!        the case kappa <= kappa0 is already covered in the public interface
    //! @param kappa history variable
    //! @return damage derivative
    double Derivative(const double kappa) const
    {
        if (kappa >= mKappaC)
            return 0.;
        double kappaScaled = GetKappaScaled(kappa);
        return -6 * mKappa0 / kappa / (mKappaC - mKappa0) * (kappaScaled * kappaScaled - kappaScaled) +
               mKappa0 / (kappa * kappa) *
                       (2 * kappaScaled * kappaScaled * kappaScaled - 3 * kappaScaled * kappaScaled + 1);
    }

private:
    double GetKappaScaled(const double kappa) const
    {
        return (kappa - mKappa0) / (mKappaC - mKappa0);
    }

    const double mKappaC;
};

} /* Constitutive */
} /* NuTo */
