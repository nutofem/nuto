#pragma once

#include "DamageLaw.h"
#include "SofteningMaterial.h"

namespace NuTo
{
namespace Constitutive
{

//! @brief damage law that exhibits no softening in the stress-strain curve and stays at peak load
//! Peerlings et al.. Gradient enhanced damage quasi-brittle materials 1996.
/*!
 *  \f[
 *  \omega = \begin{cases}
 *     0 & \text{if } \kappa < \kappa_0 \\
 *     1 - \frac{\kappa_i}{\kappa} & \text{otherwise}.
    \end{cases} \\
 *  \f]
 */
class DamageLawNoSoftening : public DamageLaw
{
public:
    DamageLawNoSoftening(double kappa0)
        : mKappa0(kappa0)
    {
    }
    DamageLawNoSoftening(Material::Softening m)
        : mKappa0(m.ft / m.E)
    {
    }

    double Damage(double kappa) const override
    {
        if (kappa < mKappa0)
            return 0.;
        return 1. - mKappa0 / kappa;
    }

    double Derivative(double kappa) const override
    {
        if (kappa < mKappa0)
            return 0.;
        return mKappa0 / (kappa * kappa);
    }

private:
    double mKappa0;
};
} /* Constitutive */
} /* NuTo */
