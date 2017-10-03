#pragma once

namespace Hygro
{


//! State equation of water from empirical research.
//! @param temperature Temperature in Kelvin.
//! @return Density of water in kg/m³
//! @remark Taken from David Jon Furbish, "Fluid Physics in Geology", 1997
double WaterDensity(const double temperature);


//! \f[ p^v = p^{vs} \mathrm{e}^{- \frac{M_w}{RT} \frac{p^c}{ρ^w}} \f]
double KelvinEquation(const double capillaryPressure, const double temperature);


//! Saturation water vapour pressure as a function of temperature.
//! @param temperature Temperature in Kelvin.
//! @return Saturation vapour pressure in MPa.
//! @remark Source: Wagner, Pruss; "International Equations for the Saturation Properties of Ordinary Water Substance."
//!         1993, DOI: [10.1063/1.555926](http://dx.doi.org/10.1063/1.555926)
double SaturationPressure(const double temperature);


const auto AirMolarMass = 28.971e-3; //! Molar mass of dry air in [kg / mol]
const auto WaterMolarMass = 18.01528e-3; // Molar mass of water [kg / mol]
} // namespace Hygro
