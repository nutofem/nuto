#include <pybind11/pybind11.h>

#include "mechanics/PDEs/HygroHelpers.h"

namespace py = pybind11;

PYBIND11_MODULE(nuto, m) {
    m.doc() = "NuTo Finite Element Library";
    m.def("DensityOfWater", &Hygro::DensityOfWater, "Density of water as a function of temperature");
}
