#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include "metamodel/NeuralNetwork.h"

namespace py = pybind11;
using namespace NuTo;

PYBIND11_MODULE(metamodel, m)
{
    m.doc() = "NuTo Metamodel python module.";

    py::class_<CallbackHandler>(m, "CallbackHandler");

    py::class_<Metamodel>(m, "Metamodel")
            .def("InitRandomNumberGenerator", &Metamodel::InitRandomNumberGenerator)
            .def("SetVerboseLevel", &Metamodel::SetVerboseLevel)
            .def("BuildTransformation", &Metamodel::BuildTransformation)
            .def("Build", &Metamodel::Build)
            .def("SolveConfidenceInterval", &Metamodel::SolveConfidenceInterval)
            .def("AppendMinMaxTransformationInput",
                 py::overload_cast<double, double>(&Metamodel::AppendMinMaxTransformationInput))
            .def("AppendMinMaxTransformationInput",
                 py::overload_cast<int, double, double>(&Metamodel::AppendMinMaxTransformationInput))
            .def("AppendMinMaxTransformationOutput",
                 py::overload_cast<double, double>(&Metamodel::AppendMinMaxTransformationOutput))
            .def("AppendMinMaxTransformationOutput",
                 py::overload_cast<int, double, double>(&Metamodel::AppendMinMaxTransformationOutput))
            .def("SetSupportPoints", &Metamodel::SetSupportPoints);

    py::class_<NeuralNetwork, Metamodel>(m, "NeuralNetwork")
            .def(py::init<std::vector<int>>())
            .def("SetBayesianTraining", &NeuralNetwork::SetBayesianTraining)
            .def("SetInitAlpha", &NeuralNetwork::SetInitAlpha)
            .def("SetAccuracyGradient", &NeuralNetwork::SetAccuracyGradient)
            .def("SetMinDeltaObjectiveBetweenRestarts", &NeuralNetwork::SetMinDeltaObjectiveBetweenRestarts)
            .def("SetMinDeltaObjectiveBayesianIteration", &NeuralNetwork::SetMinDeltaObjectiveBayesianIteration)
            .def("SetMaxFunctionCalls", &NeuralNetwork::SetMaxFunctionCalls)
            .def("SetShowSteps", &NeuralNetwork::SetShowSteps)
            .def("SetTransferFunction", &NeuralNetwork::SetTransferFunction)
            .def("SetMaxBayesianIterations", &NeuralNetwork::SetMaxBayesianIterations)
            .def("Objective", &NeuralNetwork::Objective)
            .def("Gradient", &NeuralNetwork::Gradient)
            .def("HessianFull", &NeuralNetwork::HessianFull)
            .def("GetParameters", &NeuralNetwork::GetParameters)
            .def("GetNumParameters", &NeuralNetwork::GetNumParameters)
            .def("UseDiagHessian", &NeuralNetwork::UseDiagHessian);


    py::enum_<NeuralNetwork::eTransferFunctions>(m, "eTransferFunctions")
            .value("TanSig", NeuralNetwork::eTransferFunctions::TanSig)
            .value("PureLin", NeuralNetwork::eTransferFunctions::PureLin)
            .export_values();
}
