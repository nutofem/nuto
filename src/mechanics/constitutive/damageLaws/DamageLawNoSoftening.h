#pragma once

#include <memory>
#include "mechanics/constitutive/damageLaws/DamageLaw.h"

namespace NuTo
{
namespace Constitutive
{

class DamageLawNoSoftening : public DamageLaw
{
public:
    static std::shared_ptr<DamageLaw> Create(double kappa0)
    {
        return std::shared_ptr<DamageLaw>(new DamageLawNoSoftening(kappa0));
    }

protected:
    DamageLawNoSoftening(double kappa0)
        : DamageLaw(kappa0)
    {
    }

    //! @brief protected virtual method for the damage calculation
    //!        the case kappa <= kappa0 is already covered in the public interface
    //! @param kappa history variable
    //! @return damage
    double Damage(const double kappa) const override
    {
        return 1. - mKappa0 / kappa;
    }

    //! @brief protected virtual method for the damage derivative calculation
    //!        the case kappa <= kappa0 is already covered in the public interface
    //! @param kappa history variable
    //! @return damage derivative
    double Derivative(const double kappa) const override
    {
        return mKappa0 / (kappa * kappa);
    }
};

} /* Constitutive */
} /* NuTo */
