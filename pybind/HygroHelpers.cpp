#include <pybind11/pybind11.h>

#include "mechanics/PDEs/HygroHelpers.h"
#include "mechanics/PDEs/PoreState.h"
#include "mechanics/PDEs/ConcreteMedium.h"

namespace py = pybind11;

using namespace Hygro;

PYBIND11_MODULE(hygro_helpers, m)
{
    m.doc() = "Collection of functions regarding the state of water, air and vapour.";

    m.def("DensityOfWater", &DensityOfWater, "Density of water as a function of temperature");
    m.def("SaturationPressure", &SaturationPressure, "Saturation water vapour pressure as a function of temperature");
    m.def("DynamicViscosityOfAir", &DynamicViscosityOfAir,
          "Dynamic viscosity μ_g of moist air as a function of pressures and temperature.");
    m.def("DynamicViscosityOfWater", &DynamicViscosityOfAir,
          "Dynamic viscosity μ_w of liquid water as a function temperature.");

    auto kelvinModule = m.def_submodule("KelvinEquation", "Kelvin equation");
    kelvinModule.def("VapourPressure", &KelvinEquation::VapourPressure, "Vapour pressure from the Kelvin equation.");
    kelvinModule.def("dTemperature", &KelvinEquation::dTemperature,
                     "First derivative of the Vapour pressure wrt to temperature");
    kelvinModule.def("dCapillaryPressure", &KelvinEquation::dCapillaryPressure,
                     "First derivative of the Vapour pressure wrt to capillary pressure");

    py::class_<PoreState>(m, "PoreState")
            .def(py::init<double, double, double>())
            .def_readonly("CapillaryPressure", &PoreState::CapillaryPressure)
            .def_readonly("GasPressure", &PoreState::GasPressure)
            .def_readonly("Temperature", &PoreState::Temperature)
            .def_readonly("VapourPressure", &PoreState::VapourPressure)
            .def_readonly("AirPressure", &PoreState::AirPressure)
            .def_readonly("AirDensity", &PoreState::AirDensity)
            .def_readonly("dVapourPressure_dCapillaryPressure", &PoreState::dVapourPressure_dCapillaryPressure)
            .def_readonly("VapourDensity", &PoreState::VapourDensity)
            .def_readonly("GasDensity", &PoreState::GasDensity)
            .def_readonly("GasMolarMass", &PoreState::GasMolarMass)
            .def_readonly("AirDynamicViscosity", &PoreState::AirDynamicViscosity)
            .def_readonly("WaterDensity", &PoreState::WaterDensity)
            .def_readonly("WaterDynamicViscosity", &PoreState::WaterDynamicViscosity);

    py::class_<ConcreteMedium>(m, "ConcreteMedium")
            .def(py::init<double, double, double, double>())
            .def("Porosity", &ConcreteMedium::Porosity)
            .def("Saturation", &ConcreteMedium::Saturation, "Saturation of pores with liquid water.")
            .def("DerivativeSaturation", &ConcreteMedium::DerivativeSaturation,
                 "Derivative of saturation of the pores with liquid water w.r.t. the capillary pressure.")
            .def("IntrinsicPermeability", &ConcreteMedium::IntrinsicPermeability, "Intrinsic permeability in m².")
            .def("EffectiveDiffusivity", &ConcreteMedium::EffectiveDiffusivity)
            .def("GasRelativePermeability", &ConcreteMedium::GasRelativePermeability)
            .def("WaterRelativePermeability", &ConcreteMedium::WaterRelativePermeability);
}
