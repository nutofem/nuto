#pragma once

#include "eigen3/Eigen/Core"

#include "physics/PhysicalConstantsSI.h"
#include "physics/PhysicalEquations.h"

#include "mechanics/constitutive/laws/PorousMedium.h"

#include "HygroHelpers.h"

namespace Hygro
{

namespace DryAirMassBalance
{
//! Corresponds to \f$ \mathbf{C}_{gg} \f$ in Gawin et al.
//! If the density of the air increases, naturally the mass of air in the pores increases.
Eigen::MatrixXd DensityChangeDueToGasPressure(const NuTo::PorousMedium& medium, const Eigen::MatrixXd& N,
                                              const double capillaryPressure)
{
    // input
    const double T = 273.15; // TODO: get temperature from element

    // state of dry air
    const double R = NuTo::SI::IdealGasConstant;
    const double M_a = Hygro::AirMolarMass;

    // Pores and saturation
    const double n = medium.GetPorosity();
    const double S_w = medium.Saturation(capillaryPressure);

    // weak form
    return N.transpose() * n * (1 - S_w) * M_a / (R * T) * N;
}


//! Corresponds to the first term in \f$ \mathbf{C}_{gc} \f$ in Gawin et al.
//! If the saturation increases, air is being displaced by liquid water.
Eigen::MatrixXd VariationOfSaturation(const NuTo::PorousMedium& medium, const Eigen::MatrixXd& N,
                                      const double capillaryPressure, const double gasPressure)
{
    const double temperature = 273.15; // TODO: get temperature from element

    // state of dry air
    const double R = NuTo::SI::IdealGasConstant;
    const double vapourPressure = KelvinEquation::VapourPressure(capillaryPressure, temperature);
    const double airPressure = gasPressure - vapourPressure;
    const double airDensity = airPressure * Hygro::AirMolarMass / (R * temperature);

    // derivative of degree of saturation
    const double dSw_dpc = medium.DerivativeSaturation(capillaryPressure);

    const double n = medium.GetPorosity();

    // the minus is not there in final eqs from Gawin et al., but I think it should be there
    return -N.transpose() * n * airDensity * dSw_dpc * N;
}


//! Corresponds to the second term in \f$ \mathbf{C}_{gc} \f$ in Gawin et al.
//! If the density of the air increases, naturally the mass of air in the pores increases.
Eigen::MatrixXd DensityChangeDueToCapillaryPressure(const NuTo::PorousMedium& medium, const Eigen::MatrixXd& N,
                                                    const double capillaryPressure)
{
    const double T = 273.15; // TODO: get temperature from element

    const double R = NuTo::SI::IdealGasConstant;

    const double S_w = medium.Saturation(capillaryPressure);
    const double n = medium.GetPorosity();
    const double M_a = Hygro::AirMolarMass;

    const double dpv_dpc = KelvinEquation::dCapillaryPressure(capillaryPressure, T);

    // Gawin et al. have M_w instead of M_a here, but I think it should be M_a
    return N.transpose() * n * (1 - S_w) * M_a * dpv_dpc / (R * T) * N;
}

} // namespace DryAirMassBalance

} // namespace Hygro
