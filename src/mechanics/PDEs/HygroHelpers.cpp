#include <cmath>
#include <array>

#include "HygroHelpers.h"
#include "physics/PhysicalConstantsSI.h"

using namespace Hygro;


double Hygro::WaterDensity(const double temperature)
{
    using std::pow;
    const double t = temperature - 273.15; // celsisus temperature
    const std::array<double, 6> a = {4.89e-7, -1.65e-9, 1.86e-12, 2.43e-13, -1.6e-15, 3.37e-18};
    const std::array<double, 6> b = {1.02e3, -7.74e-1, 8.77e-3, -9.21e-5, 3.35e-7, -4.4e-10};
    const double fst = b[0] + b[1] * t + b[2] * pow(t, 2) + b[3] * pow(t, 3) + b[4] * pow(t, 4) + b[5] * pow(t, 5);
    const double snd = a[0] + a[1] * t + a[2] * pow(t, 2) + a[3] * pow(t, 3) + a[4] * pow(t, 4) + a[5] * pow(t, 5);
    return fst - 1e7 * snd;
}


double Hygro::KelvinEquation(const double capillaryPressure, const double temperature)
{
    const double saturationPressure = 0.0;
    const double molarMassWater = 18.01528; // g/mol
    const double R = NuTo::SI::IdealGasConstant;
    const double waterDensity = WaterDensity(temperature);
    const double exponent = -molarMassWater * capillaryPressure / (R * temperature * waterDensity);
    return saturationPressure * std::exp(exponent);
}


double Hygro::SaturationPressure(const double temperature)
{
    const double criticalPointTemperature = 647.096; // Kelvin
    const double criticalPointPressure = 22.064; // MPa
    const double tau = 1.0 - temperature / criticalPointTemperature;
    const std::array<double, 6> a = {-7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502};
    const double polynomial = a[0] * tau + a[1] * pow(tau, 1.5) + a[2] * pow(tau, 3.0) + a[3] * pow(tau, 3.5) +
                              a[4] * pow(tau, 4.0) + a[5] * pow(tau, 7.5);
    const double exponent = criticalPointTemperature * polynomial / temperature;
    return criticalPointPressure * std::exp(exponent);
}
