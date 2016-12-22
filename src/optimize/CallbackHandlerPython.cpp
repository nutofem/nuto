// $Id$

#include <Python.h>
#include "math/FullMatrix.h"
#include "optimize/CallbackHandlerPython.h"
#include "optimize/OptimizeException.h"

namespace NuTo
{

    void  CallbackHandlerPython::SetCallbackFunctions(PyObject *args_parameters, PyObject *args_objective,PyObject *args_gradient, PyObject *args_hessian)
	{
        // check if objective routine is callable
        if (!PyCallable_Check(args_parameters)) 
		{
            throw(OptimizeException("[CallbackHandlerPython::SetCallbackFunctions] Set parameter routine not callable - check the routine"));
        }

        Py_XINCREF(args_parameters);                // Add a reference to new callback 
        Py_XDECREF(mCallbackSetParameters);         // Dispose of previous callback 
        mCallbackSetParameters = args_parameters;   // Remember new callback 

        // check of objective routine is callable
        if (!PyCallable_Check(args_objective)) 
        {
            throw(OptimizeException("[CallbackHandlerPython::SetCallbackFunctions] Objective routine not callable - check the routine"));
        }
        Py_XINCREF(args_objective);                // Add a reference to new callback 
        Py_XDECREF(mCallbackObjective);            // Dispose of previous callback 
        mCallbackObjective = args_objective;       // Remember new callback 

        if (!PyCallable_Check(args_gradient)) 
        {
            throw(OptimizeException("[CallbackHandlerPython::SetCallbackFunctions] Gradient routine not callable - check the routine"));
        }
        Py_XINCREF(args_gradient);                 // Add a reference to new callback
        Py_XDECREF(mCallbackGradient);                   // Dispose of previous callback
        mCallbackGradient = args_gradient;      // Remember new callback

        if (!PyCallable_Check(args_hessian)) 
        {
            throw(OptimizeException("[CallbackHandlerPython::SetCallbackFunctions] Hessian routine not callable - check the routine"));
        }
        Py_XINCREF(args_hessian);                 // Add a reference to new callback
        Py_XDECREF(mCallbackHessian);             // Dispose of previous callback
        mCallbackHessian = args_hessian    ;      // Remember new callback
    }

	void CallbackHandlerPython::SetParameters(const NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rParameters)
    {    
        //printf("CallbackHandlerPython::SetParameterscall SetParameters routine %p\n",mCallbackSetParameters);
		
		if (rParameters.GetNumColumns()!=1)
            throw(OptimizeException("[CallbackHandlerPython::SetParameters] Number of columns of Parameter matrix not equal to one."));

        // create python tuple as wrapper - this is not the optimal solution, better to wrap the FullMatrix directly
        PyObject *pyList, *item;

		pyList = PyList_New(rParameters.GetNumRows());
		for (int i=0; i<rParameters.GetNumRows(); i++) 
		{
			item = PyFloat_FromDouble(rParameters(i,0));
			PyList_SET_ITEM(pyList, i, item);
		}
		
		PyObject *arglist = Py_BuildValue("(O)", pyList);

        //PyObject *arglist = PyCObject_FromVoidPtr( const_cast<void*>((void*)&rParameters) , 0);
        // call external python function
        PyEval_CallObject(mCallbackSetParameters, arglist);
        Py_DECREF(arglist);
        Py_DECREF(pyList);
	}
	
	double CallbackHandlerPython::Objective()const
	{
        //printf("CallbackHandlerPython::Objective call Objective routine %p\n",mCallbackObjective);
        double objective;
        // call external python function
        PyObject *result = PyEval_CallObject(mCallbackObjective,0);

        if (result==0)
        {
            throw(OptimizeException("[CallbackHandlerPython::Objective] Objective error calling callback objective."));
        }
        objective = PyFloat_AsDouble(result);

        // decrease reference count for result
        Py_DECREF(result);
        return objective;
    }

	void CallbackHandlerPython::Gradient (NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rGradient)const
	{
	    //printf("CallbackHandler::Gradient\n");
        
        // Variante with the creation of lists, better directly wrap a full matrix
        // call external python function
        PyObject *result = PyEval_CallObject(mCallbackGradient, 0);

		if (result==0)
		{
            throw(OptimizeException("[CallbackHandlerPython::Gradient] error calling CallbackHandler::Gradient."));
		}
		
		if (!PyList_Check(result))
        {        
            throw(OptimizeException("[CallbackHandlerPython::Gradient] result is not a python list - check your python code."));
		}

        //resize gradient vector
        rGradient.Resize(PyList_GET_SIZE(result),1);

        // write the lsit back
		for (int i=0; i<rGradient.GetNumRows(); i++) 
		{
			rGradient(i,0) = PyFloat_AsDouble(PyList_GET_ITEM(result, i));
		}

        Py_DECREF(result);
    }


    void CallbackHandlerPython::Hessian(NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rHessian)const
    {
        //printf("CallbackHandler::Hessian\n");
        
        // Variante with the creation of list - better wrap directly a full matrix
        // call external python function
        PyObject *result = PyEval_CallObject(mCallbackHessian, 0);
        
        if (result==0)
        {
            throw(OptimizeException("[CallbackHandlerPython::Hessian] error calling CallbackHandler::Hessian."));
        }
        
        if (!PyList_Check(result))
        {        
            throw(OptimizeException("[CallbackHandlerPython::Hessian] result is not a python list - check your python code."));
        }
        int dim=(int)sqrt(PyList_GET_SIZE(result));
        if (dim*dim!= PyList_GET_SIZE(result))
        {
            throw(OptimizeException("[CallbackHandlerPython::Hessian] the python list for the hessian matrix should have NumParameters*NumParameters entries."));   
        }
        
        //resize Hessian matrix
        rHessian.Resize(dim,dim);
        
        // write the list back
        int entry(0);
        for (int i=0; i<dim; i++) 
        {
            for (int j=0; j<dim; j++,entry++) 
            {
                rHessian(i,j) = PyFloat_AsDouble(PyList_GET_ITEM(result, entry));
            }
        }
        
        Py_DECREF(result);
    }
    void CallbackHandlerPython::Info()const
    {
        std::cout << "CallbackHandlerPython" << std::endl;
    }
}
