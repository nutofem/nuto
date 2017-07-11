#pragma once

#include "optimize/CallbackHandler.h"

#include <Python.h>

namespace NuTo
{

//! @author Joerg F. Unger, ISM
//! @date September 2009
//! @brief Abstract class to handle callback routines
class CallbackHandlerPython : public CallbackHandler
{

public:
    CallbackHandlerPython()
        : CallbackHandler()
    {
        mCallbackSetParameters = nullptr;
        mCallbackObjective = nullptr;
        mCallbackGradient = nullptr;
        mCallbackHessian = nullptr;
    }

    void SetCallbackFunctions(PyObject* args_parameters, PyObject* args_objective, PyObject* args_gradient,
                              PyObject* args_hessian);

    void SetParameters(const Eigen::MatrixXd& rParameters) override;

    double Objective() const override;

    void Gradient(Eigen::MatrixXd& rGradient) const override;

    void Hessian(Eigen::MatrixXd& rHessian) const override;

private:
    PyObject* mCallbackSetParameters;
    PyObject* mCallbackObjective;
    PyObject* mCallbackGradient;
    PyObject* mCallbackHessian;
};
} // namespace NuTo
