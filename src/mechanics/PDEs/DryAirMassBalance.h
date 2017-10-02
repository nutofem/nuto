#pragma once

#include "eigen3/Eigen/Core"

#include "physics/PhysicalConstantsSI.h"
#include "physics/PhysicalEquations.h"

#include "mechanics/constitutive/laws/PorousMedium.h"

#include "HygroHelpers.h"

namespace Hygro
{

class DryAirMassBalance
{
public:
    //! Corresponds to \f[ \mathbf{C}_{gg} \f] in Gawin et al.
    Eigen::MatrixXd DensityChangeDueToGasPressure(const NuTo::PorousMedium& medium, const Eigen::MatrixXd& N,
                                                  const double capillaryPressure)
    {
        // input
        const double T = 298.15; // TODO: get temperature from element

        // state of dry air
        const double R = NuTo::SI::IdealGasConstant;
        const double M_a = 28.971e-3; // molar mass of dry air

        // Pores and saturation
        const double n = medium.GetPorosity();
        const double S_w = medium.Saturation(capillaryPressure);

        // weak form
        return N.transpose() * n * (1 - S_w) * M_a / (R * T) * N;
    }
};

} // namespace Hygro
