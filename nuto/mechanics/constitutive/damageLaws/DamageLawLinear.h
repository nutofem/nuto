#pragma once

#include <algorithm> // std::min
#include "nuto/mechanics/constitutive/damageLaws/DamageLaw.h"

namespace NuTo
{
namespace Constitutive
{
//! @brief linear damage law
//! Peerlings et al.. Gradient enhanced damage quasi-brittle materials 1996.
/*!
 *  \f[
 *  \omega = \begin{cases}
 *     0 & \text{if } \kappa < \kappa_0 \\
 *     \frac{\kappa_c}{\kappa} \frac{\kappa - \kappa_i}{\kappa_c - \kappa_i} & \text{otherwise}.
    \end{cases} \\
    \omega = \min(\omega, \omega_{\max})
 *  \f]
 */
class DamageLawLinear : public DamageLaw
{
public:
    DamageLawLinear(double kappa0, double kappaC, double omegaMax)
        : mKappa0(kappa0)
        , mKappaC(kappaC)
        , mOmegaMax(omegaMax)
    {
    }

    double Damage(double kappa) const override
    {
        if (kappa < mKappa0)
            return 0.;
        double damage = mKappaC / kappa * (kappa - mKappa0) / (mKappaC - mKappa0);
        return std::min(damage, mOmegaMax);
    }

    double Derivative(double kappa) const override
    {
        if (kappa < mKappa0)
            return 0.;
        if (Damage(kappa) < mOmegaMax)
            return mKappaC * mKappa0 / (kappa * kappa * (mKappaC - mKappa0));
        return 0;
    }

private:
    double mKappa0;
    double mKappaC;
    double mOmegaMax;
};
} /* Constitutive */
} /* NuTo */
