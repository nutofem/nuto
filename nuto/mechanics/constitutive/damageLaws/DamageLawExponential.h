#pragma once

#include <cmath>
#include "DamageLaw.h"
#include "SofteningMaterial.h"

namespace NuTo
{
namespace Constitutive
{
//! @brief exponential damage omega
//! Peerlings, R., De Borst, R., Brekelmans, W., Geers, M.. Gradient-enhanced damage modelling of concrete fracture.
//! Mechanics of Cohesive-frictional Materials 1998;3(4):323–342.
/*!
 *  \f[
 *  \omega = \begin{cases}
 *     0 & \text{if } \kappa < \kappa_0 \\
 *     1 - \frac{\kappa_0}{\kappa} \left(1 - \alpha + \alpha \exp \left( \frac{f_t}{g_f} (\kappa_0 - \kappa) \right)
 \right) & \text{otherwise}.
    \end{cases}
 *  \f]
 */
class DamageLawExponential : public DamageLaw
{
public:
    DamageLawExponential(double kappa0, double beta, double alpha)
        : mKappa0(kappa0)
        , mBeta(beta)
        , mAlpha(alpha)
    {
    }
    DamageLawExponential(Material::Softening m)
        : mKappa0(m.ft / m.E)
        , mBeta(m.ft / m.gf)
        , mAlpha(1 - m.fMin / m.ft)
    {
    }

    double Damage(double kappa) const override
    {
        if (kappa < mKappa0)
            return 0.;
        return 1 - mKappa0 / kappa * (1 - mAlpha + mAlpha * std::exp(mBeta * (mKappa0 - kappa)));
    }

    double Derivative(double kappa) const override
    {
        if (kappa < mKappa0)
            return 0.;
        return mKappa0 / kappa *
               ((1 / kappa + mBeta) * mAlpha * std::exp(mBeta * (mKappa0 - kappa)) + (1 - mAlpha) / kappa);
    }

private:
    double mKappa0;
    double mBeta;
    double mAlpha;
};
} /* Constitutive */
} /* NuTo */
