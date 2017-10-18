#pragma once

#include <eigen3/Eigen/Core>

#include "PoreState.h"
#include "mechanics/constitutive/laws/PorousMedium.h"

namespace Hygro
{

namespace WaterMassBalance
{

//! Corresponds to second term in \f$ \mathbf{C}_{cc} \f$ in Gawin et al.
Eigen::MatrixXd VariationOfSaturation(const PoreState& poreState, const NuTo::PorousMedium& medium,
                                      const Eigen::MatrixXd& N)
{
    const double n = medium.Porosity();
    const double waterDensity = poreState.WaterDensity;
    const double vapourDensity = poreState.VapourDensity;
    const double dSw_dpc = medium.DerivativeSaturation(poreState.CapillaryPressure);

    return N.transpose() * n * (waterDensity - vapourDensity) * dSw_dpc * N;
}

//! Corresponds to first term in \f$ \mathbf{C}_{cc} \f$ in Gawin et al.
Eigen::MatrixXd ChangeOfVapourDensity(const PoreState& poreState, const NuTo::PorousMedium& medium,
                                      const Eigen::MatrixXd& N)
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
Eigen::MatrixXd DiffusionGasPressure(const PoreState& poreState, const NuTo::PorousMedium& medium,
                                     const Eigen::MatrixXd& dN)
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
Eigen::MatrixXd DiffusionCapillaryPressure(const PoreState& poreState, const NuTo::PorousMedium& medium,
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
Eigen::MatrixXd AdvectionVapourGasPressure(const PoreState& poreState, const NuTo::PorousMedium& medium,
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
Eigen::MatrixXd AdvectionWaterGasPressure(const PoreState& poreState, const NuTo::PorousMedium& medium,
                                          const Eigen::MatrixXd& dN)
{
    const double waterDensity = poreState.WaterDensity;

    const double relPermeability = medium.WaterRelativePermeability(poreState.CapillaryPressure);
    const double intrinsicPermeability = medium.IntrinsicPermeability(poreState.CapillaryPressure);

    const double waterViscosity = poreState.WaterDynamicViscosity;

    return dN.transpose() * waterDensity * relPermeability * intrinsicPermeability * dN / waterViscosity;
}


//! Corresponds to second term in \f$ \mathbf{K}_{cc} \f$ in Gawin et al.
Eigen::MatrixXd AdvectionWaterCapillaryPressure(const PoreState& poreState, const NuTo::PorousMedium& medium,
                                                const Eigen::MatrixXd& dN)
{
    const double waterDensity = poreState.WaterDensity;

    const double relPermeability = medium.WaterRelativePermeability(poreState.CapillaryPressure);
    const double intrinsicPermeability = medium.IntrinsicPermeability(poreState.CapillaryPressure);

    const double waterViscosity = poreState.WaterDynamicViscosity;


    return -dN.transpose() * waterDensity * relPermeability * intrinsicPermeability * dN / waterViscosity;
}

} // namespace WaterMassBalance
} // namespace Hygro
