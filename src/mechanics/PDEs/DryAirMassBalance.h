#pragma once

#include "eigen3/Eigen/Core"

#include "physics/PhysicalConstantsSI.h"
#include "physics/PhysicalEquations.h"

#include "mechanics/constitutive/laws/PorousMedium.h"

#include "PoreState.h"

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
                                      const double capillaryPressure, const PoreState& poreState)
{
    const double n = medium.GetPorosity();
    const double airDensity = poreState.AirDensity;
    const double dSw_dpc = medium.DerivativeSaturation(capillaryPressure);

    // the minus is not there in final eqs from Gawin et al., but I think it should be there
    return -N.transpose() * n * airDensity * dSw_dpc * N;
}


//! Corresponds to the second term of \f$ \mathbf{C}_{gc} \f$ in Gawin et al.
//! If the density of the air increases, naturally the mass of air in the pores increases.
Eigen::MatrixXd DensityChangeDueToCapillaryPressure(const NuTo::PorousMedium& medium, const Eigen::MatrixXd& N,
                                                    const double capillaryPressure, const PoreState poreState)
{
    const double T = 273.15; // TODO: get temperature from element
    const double R = NuTo::SI::IdealGasConstant;
    const double M_a = Hygro::AirMolarMass;

    const double S_w = medium.Saturation(capillaryPressure);
    const double n = medium.GetPorosity();

    const double dpv_dpc = poreState.dVapourPressure_dCapillaryPressure;

    // Gawin et al. have M_w instead of M_a here, but I think it should be M_a
    return N.transpose() * n * (1 - S_w) * M_a * dpv_dpc / (R * T) * N;
}


//! Corresponds to first term of \f$ \mathbf{K}_{gg} \f$ in Gawin et al.
//! Contribution of advective air flow.
template <int TDim>
Eigen::MatrixXd PermeabilityMatrix(const NuTo::PorousMedium& medium, const Eigen::MatrixXd& dN,
                                   const PoreState poreState)
{
    const double airDensity = poreState.AirDensity;

    // permeability
    const double capillaryPressure = poreState.CapillaryPressure;
    const double gasPressure = poreState.GasPressure;
    const double S_w = medium.Saturation(capillaryPressure);
    const double relPermability = 1.0 - S_w;
    const Eigen::Matrix<double, TDim, TDim> intrinsicPermeability = medium.IntrinsicPermeability<TDim>(gasPressure);

    const double viscosity = poreState.AirDynamicViscosity;

    return dN.transpose() * airDensity * relPermability * intrinsicPermeability / viscosity * dN;
}


//! Corresponds to the second term of \f$ \mathbf{K}_{gg} \f$ in Gawin et al.
//! Contribution of diffusive air flow.
Eigen::MatrixXd DiffusionGasPressure(const Eigen::MatrixXd& dN, const double capillaryPressure, const double gasPressure, const PoreState poreState)
{
    const double M_a = Hygro::AirMolarMass;
    const double M_w = Hygro::WaterMolarMass;
    const double M_g = poreState.GasMolarMass;
    const double molarRatio = M_a * M_w / std::pow(M_g, 2);

    const double vapourPressure = poreState.VapourPressure;
    const double pressureRatio = vapourPressure / std::pow(gasPressure, 2);

    // TODO: diffusionTensor
    //return dN.transpose() * gasDensity * molarRatio * diffusionTensor * pressureRatio * dN;
}


//! Corresponds to \f$ \mathbf{K}_{gc} \f$ in Gawin et al.
Eigen::MatrixXd DiffusionCapillaryPressure(const Eigen::MatrixXd& dN, const double capillaryPressure, const double gasPressure, const PoreState poreState)
{
    const double M_a = Hygro::AirMolarMass;
    const double M_w = Hygro::WaterMolarMass;
    const double M_g = poreState.GasMolarMass;
    const double molarRatio = M_a * M_w / std::pow(M_g, 2);

    const double dpv_dpc = poreState.dVapourPressure_dCapillaryPressure;

    //return dN.transpose() * gasDensity * molarRatio * diffusionTensor * dpv_dpc / gasPressure * dN;
}


} // namespace DryAirMassBalance

} // namespace Hygro
