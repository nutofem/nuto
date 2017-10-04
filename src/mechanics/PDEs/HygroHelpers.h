#pragma once

namespace Hygro
{
//! State equation of water from empirical research.
//! @param temperature Temperature in Kelvin.
//! @return Density of water in kg/m³
//! @remark Taken from David Jon Furbish, "Fluid Physics in Geology", 1997
double WaterDensity(const double temperature);

namespace KelvinEquation
{
//! \f[ p^v = p^{vs} \mathrm{e}^{- \frac{M_w}{RT} \frac{p^c}{ρ^w}} \f]
//! @param capillaryPressure Capillary pressure \f$p^c\f$ in [MPa]
//! @param temperature Temperature in [K]
//! @return Vapour pressure
double VapourPressure(const double capillaryPressure, const double temperature);


//! \f[ \frac{∂p^v}{∂T} \f]
//! @param capillaryPressure Capillary pressure \f$p^c\f$ in [MPa]
//! @param temperature Temperature in [K]
//! @return Derivative wrt to temperature
double dTemperature(const double capillaryPressure, const double temperature);


//! \f[ \frac{∂p^v}{∂p^c} \f]
//! @param capillaryPressure Capillary pressure \f$p^c\f$ in [MPa]
//! @param temperature Temperature in [K]
//! @return Derivative wrt to the capillary pressure
double dCapillaryPressure(const double capillaryPressure, const double temperature);
}


//! Saturation water vapour pressure as a function of temperature.
//! @param temperature Temperature in Kelvin.
//! @return Saturation vapour pressure in MPa.
//! @remark Source: Wagner, Pruss; "International Equations for the Saturation Properties of Ordinary Water Substance."
//!         1993, DOI: [10.1063/1.555926](http://dx.doi.org/10.1063/1.555926)
double SaturationPressure(const double temperature);


//! Dynamic viscosity \f$ μ_g \f$ of moist air as a function of pressures and temperature.
//! @param airPressure Air pressure in MPa.
//! @param gasPressure Gas pressure in MPa.
//! @param temperature Temperature in Kelvin.
//! @return Dynamic viscosity in [Pa s].
//! @remark Source: Gawin et al.; "Numerical analysis of hygro-thermal behaviour and damage of concrete at
//!                 high temperature", 1999, DOI:
//! [10.1002/(SICI)1099-1484(199901)4:1<37::AID-CFM58>3.0.CO;2-S](https://dx.doi.org/10.1002/(SICI)1099-1484(199901)4:1<37::AID-CFM58>3.0.CO;2-S)
double DynamicViscosityOfAir(const double airPressure, const double gasPressure, const double temperature);


const double AirMolarMass = 28.971e-3; //!< Molar mass of dry air in [kg / mol]
const double WaterMolarMass = 18.01528e-3; //!< Molar mass of water [kg / mol]
} // namespace Hygro
