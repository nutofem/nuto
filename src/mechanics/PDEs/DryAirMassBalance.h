#pragma once

#include <eigen3/Eigen/Core>

#include "PoreState.h"
#include "PorousMedium.h"

namespace Hygro
{

//! The discretized weak forms are based on those of Gawin et al.
//! @see Gawin et al, 2006: [10.1002/nme.1615](https://dx.doi.org/10.1002/nme.1615)
//! @see Gawin et al, 2003: [10.1016/S0045-7825(03)00200-7](https://dx.doi.org/10.1016/S0045-7825(03)00200-7)
namespace DryAirMassBalance
{

//! Corresponds to \f$ \mathbf{C}_{gg} \f$ in Gawin et al.
//! If the density of the air increases, naturally the mass of air in the pores increases.
Eigen::MatrixXd DensityChangeDueToGasPressure(const PoreState& poreState, const PorousMedium& medium,
                                              const Eigen::MatrixXd& N)
{
    const double T = poreState.Temperature;

    // state of dry air
    const double R = Hygro::MolarGasConstant;
    const double M_a = Hygro::AirMolarMass;

    // Pores and saturation
    const double n = medium.Porosity();
    const double S_w = medium.Saturation(poreState.CapillaryPressure);

    // weak form
    return N.transpose() * n * (1 - S_w) * M_a / (R * T) * N;
}


//! Corresponds to the first term of \f$ \mathbf{C}_{gc} \f$ in Gawin et al.
//! If the saturation increases, air is being displaced by liquid water.
Eigen::MatrixXd VariationOfSaturation(const PoreState& poreState, const PorousMedium& medium, const Eigen::MatrixXd& N)
{
    const double n = medium.Porosity();
    const double airDensity = poreState.AirDensity;
    const double dSw_dpc = medium.DerivativeSaturation(poreState.CapillaryPressure);

    // the minus is not there in final eqs from Gawin et al., but I think it should be there
    return -N.transpose() * n * airDensity * dSw_dpc * N;
}


//! Corresponds to the second term of \f$ \mathbf{C}_{gc} \f$ in Gawin et al.
//! If the density of the air increases, naturally the mass of air in the pores increases.
Eigen::MatrixXd DensityChangeDueToCapillaryPressure(const PoreState& poreState, const PorousMedium& medium,
                                                    const Eigen::MatrixXd& N)
{
    const double T = poreState.Temperature;
    const double R = Hygro::MolarGasConstant;
    const double M_a = Hygro::AirMolarMass;

    const double n = medium.Porosity();
    const double S_w = medium.Saturation(poreState.CapillaryPressure);

    const double dpv_dpc = poreState.dVapourPressure_dCapillaryPressure;

    // Gawin et al. have M_w instead of M_a here, but I think it should be M_a
    return N.transpose() * n * (1 - S_w) * M_a * dpv_dpc / (R * T) * N;
}


//! Corresponds to first term of \f$ \mathbf{K}_{gg} \f$ in Gawin et al.
//! Contribution of advective air flow.
Eigen::MatrixXd AdvectionGasPressure(const PoreState& poreState, const PorousMedium& medium, const Eigen::MatrixXd& dN)
{
    const double airDensity = poreState.AirDensity;

    const double relPermeability = medium.GasRelativePermeability(poreState.CapillaryPressure);

    const double gasPressure = poreState.GasPressure;
    const double intrinsicPermeability = medium.IntrinsicPermeability(gasPressure);

    const double viscosity = poreState.AirDynamicViscosity;

    return dN.transpose() * airDensity * relPermeability * intrinsicPermeability / viscosity * dN;
}


//! Corresponds to the second term of \f$ \mathbf{K}_{gg} \f$ in Gawin et al.
//! Contribution of diffusive air flow.
Eigen::MatrixXd DiffusionGasPressure(const PoreState& poreState, const PorousMedium& medium, const Eigen::MatrixXd& dN)
{
    const double gasDensity = poreState.GasDensity;

    const double M_a = Hygro::AirMolarMass;
    const double M_w = Hygro::WaterMolarMass;
    const double M_g = poreState.GasMolarMass;
    const double molarRatio = M_a * M_w / std::pow(M_g, 2);

    const double gasPressure = poreState.GasPressure;
    const double vapourPressure = poreState.VapourPressure;
    const double pressureRatio = vapourPressure / std::pow(gasPressure, 2);

    const double capillaryPressure = poreState.CapillaryPressure;
    const double temperature = poreState.Temperature;
    const double diffusionTensor = medium.EffectiveDiffusivity(capillaryPressure, gasPressure, temperature);

    return dN.transpose() * gasDensity * molarRatio * diffusionTensor * pressureRatio * dN;
}


//! Corresponds to \f$ \mathbf{K}_{gc} \f$ in Gawin et al.
Eigen::MatrixXd DiffusionCapillaryPressure(const PoreState& poreState, const PorousMedium& medium,
                                           const Eigen::MatrixXd& dN)
{
    const double gasDensity = poreState.GasDensity;

    const double M_a = Hygro::AirMolarMass;
    const double M_w = Hygro::WaterMolarMass;
    const double M_g = poreState.GasMolarMass;
    const double molarRatio = M_a * M_w / std::pow(M_g, 2);

    const double temperature = poreState.Temperature;
    const double capillaryPressure = poreState.CapillaryPressure;
    const double gasPressure = poreState.GasPressure;
    const double diffusionTensor = medium.EffectiveDiffusivity(capillaryPressure, gasPressure, temperature);

    const double dpv_dpc = poreState.dVapourPressure_dCapillaryPressure;

    return dN.transpose() * gasDensity * molarRatio * diffusionTensor * dpv_dpc / gasPressure * dN;
}

//! Corresponds to \f$ \mathbf{f}_g \f$ in Gawin et al.
Eigen::VectorXd AdvectiveGravityLoad(const PoreState& poreState, const PorousMedium& medium, const Eigen::MatrixXd& dN)
{
    const double airDensity = poreState.AirDensity;

    const double relPermeability = medium.GasRelativePermeability(poreState.CapillaryPressure);

    const double gasPressure = poreState.GasPressure;
    const double intrinsicPermeability = medium.IntrinsicPermeability(gasPressure);

    const double viscosity = poreState.AirDynamicViscosity;

    const double gasDensity = poreState.GasDensity;

    const int dimension = dN.rows();
    Eigen::VectorXd g = Eigen::VectorXd::Zero(dimension);
    g.tail(1) = 9.80665 * Eigen::VectorXd::Ones(1);

    return dN.transpose() * airDensity * relPermeability * intrinsicPermeability / viscosity * gasDensity * g;
}

} // namespace DryAirMassBalance

} // namespace Hygro
