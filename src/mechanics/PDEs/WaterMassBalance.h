#pragma once

#include <eigen3/Eigen/Core>

#include "PoreState.h"
#include "PorousMedium.h"

namespace Hygro
{

//! The discretized weak forms are based on those of Gawin et al.
//! @see Gawin et al, 2006: [10.1002/nme.1615](https://dx.doi.org/10.1002/nme.1615)
//! @see Gawin et al, 2003: [10.1016/S0045-7825(03)00200-7](https://dx.doi.org/10.1016/S0045-7825(03)00200-7)
namespace WaterMassBalance
{

//! Corresponds to second term in \f$ \mathbf{C}_{cc} \f$ in Gawin et al.
Eigen::MatrixXd VariationOfSaturation(const PoreState& poreState, const PorousMedium& medium, const Eigen::MatrixXd& N)
{
    const double n = medium.Porosity();
    const double waterDensity = poreState.WaterDensity;
    const double vapourDensity = poreState.VapourDensity;
    const double dSw_dpc = medium.DerivativeSaturation(poreState.CapillaryPressure);

    return N.transpose() * n * (waterDensity - vapourDensity) * dSw_dpc * N;
}

//! Corresponds to first term in \f$ \mathbf{C}_{cc} \f$ in Gawin et al.
Eigen::MatrixXd ChangeOfVapourDensity(const PoreState& poreState, const PorousMedium& medium, const Eigen::MatrixXd& N)
{
    const double n = medium.Porosity();
    const double S_w = medium.Saturation(poreState.CapillaryPressure);
    const double M_w = Hygro::WaterMolarMass;
    const double dpv_dpc = poreState.dVapourPressure_dCapillaryPressure;
    const double R = Hygro::MolarGasConstant;
    const double T = poreState.Temperature;

    return N.transpose() * n * (1 - S_w) * M_w * dpv_dpc / (R * T) * N;
}

//! Corresponds to first term in \f$ \mathbf{K}_{cg} \f$ in Gawin et al.
Eigen::MatrixXd DiffusionGasPressure(const PoreState& poreState, const PorousMedium& medium, const Eigen::MatrixXd& dN)
{
    const double airDensity = poreState.AirDensity;

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

    return -dN.transpose() * airDensity * molarRatio * diffusionTensor * pressureRatio * dN;
}

//! Corresponds to first term in \f$ \mathbf{K}_{cc} \f$ in Gawin et al.
Eigen::MatrixXd DiffusionCapillaryPressure(const PoreState& poreState, const PorousMedium& medium,
                                           const Eigen::MatrixXd& dN)
{
    const double gasDensity = poreState.GasDensity;

    const double M_a = Hygro::AirMolarMass;
    const double M_w = Hygro::WaterMolarMass;
    const double M_g = poreState.GasMolarMass;
    const double molarRatio = M_a * M_w / std::pow(M_g, 2);

    const double gasPressure = poreState.GasPressure;
    const double capillaryPressure = poreState.CapillaryPressure;
    const double temperature = poreState.Temperature;
    const double diffusionTensor = medium.EffectiveDiffusivity(capillaryPressure, gasPressure, temperature);

    const double dpv_dpc = poreState.dVapourPressure_dCapillaryPressure;

    return dN.transpose() * gasDensity * molarRatio * diffusionTensor * dpv_dpc * dN / gasPressure;
}

//! Corresponds to second term in \f$ \mathbf{K}_{cg} \f$ in Gawin et al.
Eigen::MatrixXd AdvectionVapourGasPressure(const PoreState& poreState, const PorousMedium& medium,
                                           const Eigen::MatrixXd& dN)
{
    const double vapourDensity = poreState.VapourDensity;

    const double gasPressure = poreState.GasPressure;
    const double intrinsicPermeability = medium.IntrinsicPermeability(gasPressure);
    const double relPermability = medium.GasRelativePermeability(poreState.CapillaryPressure);

    const double gasViscosity = poreState.AirDynamicViscosity;
    return dN.transpose() * vapourDensity * relPermability * intrinsicPermeability * dN / gasViscosity;
}


//! Corresponds to third term in \f$ \mathbf{K}_{cg} \f$ in Gawin et al.
Eigen::MatrixXd AdvectionWaterGasPressure(const PoreState& poreState, const PorousMedium& medium,
                                          const Eigen::MatrixXd& dN)
{
    const double waterDensity = poreState.WaterDensity;

    const double relPermeability = medium.WaterRelativePermeability(poreState.CapillaryPressure);
    const double intrinsicPermeability = medium.IntrinsicPermeability(poreState.CapillaryPressure);

    const double waterViscosity = poreState.WaterDynamicViscosity;

    return dN.transpose() * waterDensity * relPermeability * intrinsicPermeability * dN / waterViscosity;
}


//! Corresponds to second term in \f$ \mathbf{K}_{cc} \f$ in Gawin et al.
Eigen::MatrixXd AdvectionWaterCapillaryPressure(const PoreState& poreState, const PorousMedium& medium,
                                                const Eigen::MatrixXd& dN)
{
    const double waterDensity = poreState.WaterDensity;

    const double relPermeability = medium.WaterRelativePermeability(poreState.CapillaryPressure);
    const double intrinsicPermeability = medium.IntrinsicPermeability(poreState.CapillaryPressure);

    const double waterViscosity = poreState.WaterDynamicViscosity;


    return -dN.transpose() * waterDensity * relPermeability * intrinsicPermeability * dN / waterViscosity;
}


//! Corresponds to \f$ \mathbf{f}_c \f$ in Gawin et al.
Eigen::VectorXd AdvectiveGravityLoad(const PoreState& poreState, const PorousMedium& medium, const Eigen::MatrixXd& dN)
{
    const double vapourDensity = poreState.VapourDensity;
    const double waterDensity = poreState.WaterDensity;
    const double gasDensity = poreState.GasDensity;

    const double relGasPermeability = medium.GasRelativePermeability(poreState.CapillaryPressure);
    const double relWaterPermeability = medium.WaterRelativePermeability(poreState.CapillaryPressure);

    const double gasPressure = poreState.GasPressure;
    const double intrinsicPermeability = medium.IntrinsicPermeability(gasPressure);

    const double gasViscosity = poreState.AirDynamicViscosity;
    const double waterViscosity = poreState.WaterDynamicViscosity;

    const int dimension = dN.rows();
    Eigen::VectorXd g = Eigen::VectorXd::Zero(dimension);
    g.tail(1) = 9.80665 * Eigen::VectorXd::Ones(1);

    const auto vapourTerm = vapourDensity * relGasPermeability * intrinsicPermeability / gasViscosity * gasDensity * g;
    const auto liquidTerm =
            waterDensity * relWaterPermeability * intrinsicPermeability / waterViscosity * waterDensity * g;
    return dN.transpose() * (vapourTerm + liquidTerm);
    // TOOD: boundary terms
}

} // namespace WaterMassBalance
} // namespace Hygro
