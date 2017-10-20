#include <cmath>
#include <array>
#include <tuple>

#include "HygroHelpers.h"

using namespace Hygro;

double Hygro::DensityOfWater(const double temperature)
{
    using std::pow;
    const double t = temperature - 273.15; // celsisus temperature
    const std::array<double, 6> a = {4.89e-7, -1.65e-9, 1.86e-12, 2.43e-13, -1.6e-15, 3.37e-18};
    const std::array<double, 6> b = {1.02e3, -7.74e-1, 8.77e-3, -9.21e-5, 3.35e-7, -4.4e-10};
    const double fst = b[0] + b[1] * t + b[2] * pow(t, 2) + b[3] * pow(t, 3) + b[4] * pow(t, 4) + b[5] * pow(t, 5);
    const double snd = a[0] + a[1] * t + a[2] * pow(t, 2) + a[3] * pow(t, 3) + a[4] * pow(t, 4) + a[5] * pow(t, 5);
    return fst - 1e7 * snd;
}


std::tuple<double, double, double> KelvinEq(const double capillaryPressure, const double temperature)
{
    const double R = MolarGasConstant;
    const double M_w = WaterMolarMass;
    const double waterDensity = DensityOfWater(temperature);
    const double factor = M_w / (R * waterDensity);

    const double T = temperature;
    const double saturationPressure= SaturationPressure(temperature);

    const double exponent = -factor * capillaryPressure / T;
    const double p_v = saturationPressure * std::exp(exponent);
    const double dpv_dpc = -factor * saturationPressure * std::exp(exponent) / T;
    const double dpv_dt = factor * saturationPressure * capillaryPressure * std::exp(exponent) / (T * T);

    return {p_v, dpv_dpc, dpv_dt};
}


double Hygro::KelvinEquation::VapourPressure(const double capillaryPressure, const double temperature)
{
    return std::get<0>(KelvinEq(capillaryPressure, temperature));
}


double Hygro::KelvinEquation::dCapillaryPressure(const double capillaryPressure, const double temperature)
{
    return std::get<1>(KelvinEq(capillaryPressure, temperature));
}


double Hygro::KelvinEquation::dTemperature(const double capillaryPressure, const double temperature)
{
    return std::get<2>(KelvinEq(capillaryPressure, temperature));
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


double Hygro::DynamicViscosityOfAir(const double airPressure, const double gasPressure, const double temperature)
{
    const double T0 = 273.15;

    const double alpha_v = 3.53e-8;
    const double mu_v0 = 8.85e-6;
    const double mu_v = mu_v0 + alpha_v * (temperature - T0);

    const double mu_a0 = 17.17e-6;
    const double alpha_a = 4.73e-8;
    const double beta_a = 2.22e-11;
    const double mu_a = mu_a0 + alpha_a * (temperature - T0) + beta_a * std::pow(temperature - T0, 2);

    const double moleFraction = airPressure / gasPressure;
    return mu_v + (mu_a - mu_v) * std::pow(moleFraction, 0.608);
}

double Hygro::DynamicViscosityOfWater(const double temperature)
{
    const double T_offset = 229.0;
    const double exponent = -1.562;
    const double factor = 0.6612;
    return factor * std::pow(temperature - T_offset, exponent);
}
