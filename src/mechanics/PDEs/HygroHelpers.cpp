#include <cmath>
#include <array>

#include "HygroHelpers.h"
#include "physics/PhysicalConstantsSI.h"
#include <boost/units/systems/si/codata/typedefs.hpp>
#include <boost/units/systems/si/codata/physico-chemical_constants.hpp>
#include <boost/units/systems/si/prefixes.hpp>
#include <boost/units/systems/si/mass_density.hpp>
#include <boost/units/systems/si/dimensionless.hpp>
#include <boost/units/systems/si/pressure.hpp>

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

using namespace boost::units;
using Pressure = quantity<si::pressure>;
using Density = quantity<si::mass_density>;
using Temperature = quantity<si::temperature>;
using Dimensionless = quantity<si::dimensionless>;
using PressurePerTemp = divide_typeof_helper<Pressure, Temperature>::type;

std::tuple<Pressure, Dimensionless, PressurePerTemp> KelvinEq(const double capillaryPressure,
                                                                    const double temperature)
{
    const auto MPa = si::mega * si::pascal;

    const auto R = si::constants::codata::R;
    const auto M_w = WaterMolarMass * si::kilogram / si::mole;
    const Density waterDensity(WaterDensity(temperature) * si::kilogram_per_cubic_meter);
    const auto factor = M_w / (R * waterDensity);

    const Temperature T = temperature * si::kelvin;
    const Pressure saturationPressure(SaturationPressure(temperature) * MPa);
    const Pressure p_c(capillaryPressure * MPa);

    const Dimensionless exponent = -factor * p_c / T;
    const Pressure p_v = saturationPressure * std::exp(exponent);
    const Dimensionless dpv_dpc = -factor * saturationPressure * std::exp(exponent) / T;
    const PressurePerTemp dpv_dt = factor * saturationPressure * p_c * std::exp(exponent) / (T * T);

    return {p_v, dpv_dpc, dpv_dt};
}


double Hygro::KelvinEquation::VapourPressure(const double capillaryPressure, const double temperature)
{
    return std::get<0>(KelvinEq(capillaryPressure, temperature)).value() * 1e-6; // MPa
}


double Hygro::KelvinEquation::dCapillaryPressure(const double capillaryPressure, const double temperature)
{
    return std::get<1>(KelvinEq(capillaryPressure, temperature));
}


double Hygro::KelvinEquation::dTemperature(const double capillaryPressure, const double temperature)
{
    return std::get<2>(KelvinEq(capillaryPressure, temperature)).value() * 1e-6; // MPa
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
