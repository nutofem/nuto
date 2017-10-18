#pragma once

#include "HygroHelpers.h"

namespace Hygro
{

class PoreState
{
public:
    PoreState(const double capillaryPressure, const double gasPressure, const double temperature)
        : CapillaryPressure(capillaryPressure)
        , GasPressure(gasPressure)
        , Temperature(temperature)
        , VapourPressure(KelvinEquation::VapourPressure(capillaryPressure, temperature))
        , AirPressure(gasPressure - VapourPressure)
        , AirDensity(AirPressure * Hygro::AirMolarMass / (R * temperature))
        , dVapourPressure_dCapillaryPressure(KelvinEquation::dCapillaryPressure(capillaryPressure, temperature))
        , VapourDensity(VapourPressure * Hygro::WaterMolarMass / (R * temperature))
        , GasDensity(AirDensity + VapourDensity)
        , GasMolarMass(1.0 / (VapourDensity / (GasDensity * Hygro::WaterMolarMass) +
                              AirDensity / (GasDensity * Hygro::AirMolarMass)))
        , AirDynamicViscosity(DynamicViscosityOfAir(AirPressure, GasPressure, Temperature))
        , WaterDensity(DensityOfWater(Temperature))
        , WaterDynamicViscosity(DynamicViscosityOfWater(Temperature))
    {
    }

    const double CapillaryPressure;
    const double GasPressure;
    const double Temperature;
    const double VapourPressure;
    const double AirPressure;
    const double AirDensity;
    const double dVapourPressure_dCapillaryPressure;
    const double VapourDensity;
    const double GasDensity;
    const double GasMolarMass;
    const double AirDynamicViscosity;
    const double WaterDensity;
    const double WaterDynamicViscosity;

private:
    static constexpr double R = Hygro::MolarGasConstant;
};

} // namespace Hygro
