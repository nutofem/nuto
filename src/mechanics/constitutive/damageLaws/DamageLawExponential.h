#pragma once

#include <memory>
#include <cmath>
#include "mechanics/constitutive/damageLaws/DamageLaw.h"

namespace NuTo
{
namespace Constitutive
{

class DamageLawExponential : public DamageLaw
{
public:
    //! @brief Create damage law with exponential softening
    //! @param kappa0 initial kappa
    //! @param beta represents tensile strength / fracture energy
    //! @param alpha (optional) max damage
    static std::shared_ptr<DamageLaw> Create(double kappa0, double beta, double alpha = 1)
    {
        return std::shared_ptr<DamageLaw>(new DamageLawExponential(kappa0, beta, alpha));
    }

protected:
    DamageLawExponential(double kappa0, double beta, double alpha)
        : DamageLaw(kappa0)
        , mBeta(beta)
        , mAlpha(alpha)
    {
    }

    //! @brief protected virtual method for the damage calculation
    //!        the case kappa <= kappa0 is already covered in the public interface
    //! @param kappa history variable
    //! @return damage
    double Damage(const double kappa) const
    {
        return 1 - mKappa0 / kappa * (1 - mAlpha + mAlpha * std::exp(mBeta * (mKappa0 - kappa)));
    }

    //! @brief protected virtual method for the damage derivative calculation
    //!        the case kappa <= kappa0 is already covered in the public interface
    //! @param kappa history variable
    //! @return damage derivative
    double Derivative(const double kappa) const
    {

        return mKappa0 / kappa *
               ((1 / kappa + mBeta) * mAlpha * std::exp(mBeta * (mKappa0 - kappa)) + (1 - mAlpha) / kappa);
    }

private:
    const double mBeta;
    const double mAlpha;
};

} /* Constitutive */
} /* NuTo */
