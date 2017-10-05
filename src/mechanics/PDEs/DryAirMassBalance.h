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


//! Corresponds to the first term of \f$ \mathbf{C}_{gc} \f$ in Gawin et al.
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


//! Corresponds to the second term of \f$ \mathbf{C}_{gc} \f$ in Gawin et al.
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


//! Corresponds to first term of \f$ \mathbf{K}_{gg} \f$ in Gawin et al.
//! Contribution of advective air flow.
template <int TDim>
Eigen::MatrixXd PermeabilityMatrix(const NuTo::PorousMedium& medium, const Eigen::MatrixXd& dN,
                                   const double capillaryPressure, const double gasPressure)
{
    const double temperature = 273.15; // TODO: get temperature from element

    // state of dry air
    const double R = NuTo::SI::IdealGasConstant;
    const double vapourPressure = KelvinEquation::VapourPressure(capillaryPressure, temperature);
    const double airPressure = gasPressure - vapourPressure;
    const double airDensity = airPressure * Hygro::AirMolarMass / (R * temperature);

    // permeability
    const double S_w = medium.Saturation(capillaryPressure);
    const double relPermability = 1.0 - S_w;
    const Eigen::Matrix<double, TDim, TDim> intrinsicPermeability = medium.IntrinsicPermeability<TDim>(gasPressure);

    const double viscosity = Hygro::DynamicViscosityOfAir(airPressure, gasPressure, temperature);

    return dN.transpose() * airDensity * relPermability * intrinsicPermeability / viscosity * dN;
}

Eigen::MatrixXd DiffusionMatrix(const Eigen::MatrixXd& dN, const double capillaryPressure, const double gasPressure)
{
    const double temperature = 273.15; // TODO: get temperature from element
    const double R = NuTo::SI::IdealGasConstant;
    const double M_a = Hygro::AirMolarMass;
    const double M_w = Hygro::WaterMolarMass;

    const double vapourPressure = KelvinEquation::VapourPressure(capillaryPressure, temperature);
    const double vapourDensity = vapourPressure * M_w / (R * temperature);
    const double airPressure = gasPressure - vapourPressure;
    const double airDensity = airPressure * Hygro::AirMolarMass / (R * temperature);
    const double gasDensity = airDensity + vapourDensity;
    const double M_g = 1.0 / (vapourDensity / (gasDensity * M_w) + airDensity / (gasDensity * M_a));

    const double molarRatio = M_a * M_w / std::pow(M_g, 2);
    const double pressureRatio = vapourPressure / std::pow(gasPressure, 2);

    // TODO: diffusionTensor
    //return dN.transpose() * gasDensity * molarRatio * diffusionTensor * pressureRatio * dN;
}

} // namespace DryAirMassBalance

} // namespace Hygro
