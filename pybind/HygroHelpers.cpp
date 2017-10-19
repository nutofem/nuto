#include <pybind11/pybind11.h>

#include "mechanics/PDEs/HygroHelpers.h"

namespace py = pybind11;

PYBIND11_MODULE(hygro_helpers, m)
{
    m.doc() = "Collection of functions regarding the state of water, air and vapour.";

    m.def("DensityOfWater", &Hygro::DensityOfWater, "Density of water as a function of temperature");
    m.def("SaturationPressure", &Hygro::SaturationPressure,
          "Saturation water vapour pressure as a function of temperature");
    m.def("DynamicViscosityOfAir", &Hygro::DynamicViscosityOfAir,
          "Dynamic viscosity μ_g of moist air as a function of pressures and temperature.");
    m.def("DynamicViscosityOfWater", &Hygro::DynamicViscosityOfAir,
          "Dynamic viscosity μ_w of liquid water as a function temperature.");

    auto kelvinModule = m.def_submodule("KelvinEquation", "Kelvin equation");
    kelvinModule.def("VapourPressure", &Hygro::KelvinEquation::VapourPressure,
                     "Vapour pressure from the Kelvin equation.");
    kelvinModule.def("dTemperature", &Hygro::KelvinEquation::dTemperature,
                     "First derivative of the Vapour pressure wrt to temperature");
    kelvinModule.def("dCapillaryPressure", &Hygro::KelvinEquation::dCapillaryPressure,
                     "First derivative of the Vapour pressure wrt to capillary pressure");
}
