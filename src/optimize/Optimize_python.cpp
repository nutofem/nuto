// $Id$

#include "optimize/Optimizer.h"
namespace NuTo
{
/*
    void Optimizer::set_callbacks(PyObject *args_objective,PyObject *args_gradient)
    {
        // check of objective routine is callable
        if (!PyCallable_Check(args_objective)) 
        {
            exit(0);
        }
        Py_XINCREF(args_objective);                // Add a reference to new callback 
        Py_XDECREF(callback_objective);            // Dispose of previous callback 
        callback_objective = args_objective;       // Remember new callback 

        if (!PyCallable_Check(args_gradient)) 
        {
            exit(0);
        }
        Py_XINCREF(args_gradient);                 // Add a reference to new callback
        Py_XDECREF(callback_gradient);                   // Dispose of previous callback
        callback_gradient = args_gradient;      // Remember new callback
    }

    double Optimizer::call_callback_objective_from_C()
    {
        double objective;
        // create python object as wrapper from std:vector<double>
        PyObject *arglist = Py_BuildValue("(dd)",Parameters[0],Parameters[1]);
        
        // call external python function
        printf("object routine %p\n",callback_objective);
        PyObject *result = PyObject_CallObject(callback_objective, arglist);
        if (result==0)
        {
            printf("error calling callback objective");
            exit(0);
        }
        objective = PyFloat_AsDouble(result);

        // decrease reference count for arglist
        Py_DECREF(arglist);
        // decrease reference count for result
        Py_DECREF(result);
        return objective;
    }

    void Optimizer::call_callback_gradient_from_C(std::vector<double>& gradient)
    {
        // create python object as wrapper from std:vector<double>
        PyObject *pygradient; //= wrap_to_python(gradient) ??
        
        // build argument list for call of python function
        PyObject *arglist = Py_BuildValue("(o)", pygradient);
        
        // call external python function
        PyObject *result = PyObject_CallObject(callback_gradient, arglist);
        if (result==0)
        {
            printf("error calling callback gradient");
            exit(0);
        }
        
        //no result required, since it is written in the vector externally given
        
        // decrease reference count for arglist
        Py_DECREF(arglist);
        // decrease reference count for result
        Py_DECREF(result);
    }
    
    void Optimizer::optimize()
    {
        double lambda = 1.0;
        for (int theCycle=0; theCycle<10;theCycle++ )
        {
            // calculate objective in python routine
            Objective = call_callback_objective_from_C();

            // calculate gradient of objective function in python routine
            call_callback_gradient_from_C(Gradient);
            
            // modify parameters in negative gradient direction with a certain step length
            for (unsigned int count=0; count<Parameters.size(); count++)
                Parameters[count]-=lambda*Gradient[count];
        }
    
    }
*/    
    }
