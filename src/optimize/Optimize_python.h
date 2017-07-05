#pragma once

#include <Python.h>


namespace NuTo
{

class Optimizer
{
public:
    Optimizer(unsigned int numParameters)
    {
        /*        callback_objective=0;
                callback_gradient=0;
                Parameters.resize(numParameters);
                Gradient.resize(numParameters);
        */
    }
    /*
        void set_callbacks(PyObject *args_objective,PyObject *args_gradient);
        double call_callback_objective_from_C();
        void call_callback_gradient_from_C(std::vector<double>& gradient);
        void optimize();
    */

private:
    /*
        PyObject *callback_objective;
        PyObject *callback_gradient;
        std::vector<double> Parameters;
        double Objective;
        std::vector<double> Gradient;
    */
};
} // namespace NuTo
