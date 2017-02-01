#include <Python.h>
#include "optimize/CallbackHandlerPython.h"
#include "optimize/OptimizeException.h"

using namespace NuTo;

void CallbackHandlerPython::SetCallbackFunctions(
        PyObject* args_parameters, PyObject* args_objective, PyObject* args_gradient, PyObject* args_hessian)
{
    // check if mObjective routine is callable
    if (!PyCallable_Check(args_parameters))
    {
        throw(OptimizeException(__PRETTY_FUNCTION__, "Set parameter routine not callable - check the routine"));
    }

    Py_XINCREF(args_parameters); // Add a reference to new callback
    Py_XDECREF(mCallbackSetParameters); // Dispose of previous callback
    mCallbackSetParameters = args_parameters; // Remember new callback

    // check of mObjective routine is callable
    if (!PyCallable_Check(args_objective))
    {
        throw(OptimizeException(__PRETTY_FUNCTION__, "Objective routine not callable - check the routine"));
    }
    Py_XINCREF(args_objective); // Add a reference to new callback
    Py_XDECREF(mCallbackObjective); // Dispose of previous callback
    mCallbackObjective = args_objective; // Remember new callback

    if (!PyCallable_Check(args_gradient))
    {
        throw(OptimizeException(__PRETTY_FUNCTION__, "Gradient routine not callable - check the routine"));
    }
    Py_XINCREF(args_gradient); // Add a reference to new callback
    Py_XDECREF(mCallbackGradient); // Dispose of previous callback
    mCallbackGradient = args_gradient; // Remember new callback

    if (!PyCallable_Check(args_hessian))
    {
        throw(OptimizeException(__PRETTY_FUNCTION__, "Hessian routine not callable - check the routine"));
    }
    Py_XINCREF(args_hessian); // Add a reference to new callback
    Py_XDECREF(mCallbackHessian); // Dispose of previous callback
    mCallbackHessian = args_hessian; // Remember new callback
}


void CallbackHandlerPython::SetParameters(const Eigen::MatrixXd& rParameters)
{
    if (rParameters.cols() != 1)
        throw(OptimizeException(__PRETTY_FUNCTION__, "Number of columns of Parameter matrix not equal to one."));

    // create python tuple as wrapper - this is not the optimal solution, better to wrap the Eigen matrix directly
    PyObject *pyList, *item;

    pyList = PyList_New(rParameters.rows());
    for (int i = 0; i < rParameters.rows(); i++)
    {
        item = PyFloat_FromDouble(rParameters(i, 0));
        PyList_SET_ITEM(pyList, i, item);
    }

    PyObject* arglist = Py_BuildValue("(O)", pyList);

    // call external python function
    PyEval_CallObject(mCallbackSetParameters, arglist);
    Py_DECREF(arglist);
    Py_DECREF(pyList);
}


double CallbackHandlerPython::Objective() const
{
    double objective;
    // call external python function
    PyObject* result = PyEval_CallObject(mCallbackObjective, 0);

    if (result == 0)
    {
        throw(OptimizeException(__PRETTY_FUNCTION__, "Objective error calling callback mObjective."));
    }
    objective = PyFloat_AsDouble(result);

    // decrease reference count for result
    Py_DECREF(result);
    return objective;
}


void CallbackHandlerPython::Gradient(Eigen::MatrixXd& rGradient) const
{
    // Variant with the creation of lists, better directly wrap a full matrix
    // call external python function
    PyObject* result = PyEval_CallObject(mCallbackGradient, 0);

    if (result == 0)
    {
        throw(OptimizeException(__PRETTY_FUNCTION__, "Error calling CallbackHandler::Gradient."));
    }

    if (!PyList_Check(result))
    {
        throw(OptimizeException(__PRETTY_FUNCTION__, "Result is not a python list - check your python code."));
    }

    rGradient.resize(PyList_GET_SIZE(result), 1);

    // write the lsit back
    for (int i = 0; i < rGradient.rows(); i++)
    {
        rGradient(i, 0) = PyFloat_AsDouble(PyList_GET_ITEM(result, i));
    }

    Py_DECREF(result);
}


void CallbackHandlerPython::Hessian(Eigen::MatrixXd& rHessian) const
{
    // Variant with the creation of list - better wrap directly a full matrix
    // call external python function
    PyObject* result = PyEval_CallObject(mCallbackHessian, 0);

    if (result == 0)
    {
        throw(OptimizeException(__PRETTY_FUNCTION__, "Error calling CallbackHandler::Hessian."));
    }

    if (!PyList_Check(result))
    {
        throw(OptimizeException(__PRETTY_FUNCTION__, "Result is not a python list - check your python code."));
    }
    int dim = (int)sqrt(PyList_GET_SIZE(result));
    if (dim * dim != PyList_GET_SIZE(result))
    {
        throw(OptimizeException(__PRETTY_FUNCTION__,
                "The python list for the hessian matrix should have NumParameters*NumParameters entries."));
    }

    rHessian.resize(dim, dim);

    // write the list back
    int entry(0);
    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j < dim; j++, entry++)
        {
            rHessian(i, j) = PyFloat_AsDouble(PyList_GET_ITEM(result, entry));
        }
    }

    Py_DECREF(result);
}


void CallbackHandlerPython::Info() const { std::cout << "CallbackHandlerPython" << std::endl; }
