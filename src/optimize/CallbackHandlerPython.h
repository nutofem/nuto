#pragma once

// parent
#include "optimize/CallbackHandler.h"

#include <Python.h>

namespace NuTo
{

//! @author Joerg F. Unger, ISM
//! @date September 2009
//! @brief ... abstract class to handle callback routines
class CallbackHandlerPython : public CallbackHandler 
{
#ifdef ENABLE_SERIALIZATION
	friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
    CallbackHandlerPython() : CallbackHandler()
    {
	    mCallbackSetParameters=0;
	    mCallbackObjective=0;
	    mCallbackGradient=0;
        mCallbackHessian=0;
    }
	
#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
    	ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(CallbackHandler);
    	mCallbackSetParameters = 0;
    	mCallbackObjective = 0;
    	mCallbackGradient = 0;
    	mCallbackHessian = 0;
    	/* the serializaion of the function pointers in python makes no sense
    	 * since in the restoring phase in another pointer session the python function might
    	 * be located at a totally different position

    	   & BOOST_SERIALIZATION_NVP(mCallbackSetParameters)
           & BOOST_SERIALIZATION_NVP(mCallbackObjective)
    	   & BOOST_SERIALIZATION_NVP(mCallbackGradient)
    	   & BOOST_SERIALIZATION_NVP(mCallbackHessian);
    	 */
    }
#endif // ENABLE_SERIALIZATION

    void  SetCallbackFunctions(PyObject *args_parameters, PyObject *args_objective, PyObject *args_gradient, PyObject *args_hessian);

	void SetParameters(const NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rParameters);
	
	double Objective()const;

	void Gradient (NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rGradient)const;
	
    void Hessian(NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rHessian)const;

    void Info()const;

#ifdef ENABLE_SERIALIZATION
    //! @brief ... restore the object from a file
    //! @param filename ... filename
    //! @param aType ... type of file, either BINARY, XML or TEXT
    //! @brief ... save the object to a file
    virtual void Restore (const std::string &filename, std::string rType )
    {

    }

	//  @brief this routine has to be implemented in the final derived classes, which are no longer abstract
    //! @param filename ... filename
    //! @param aType ... type of file, either BINARY, XML or TEXT
	virtual void Save (const std::string &filename, std::string rType )const
	{

	}
#endif // ENABLE_SERIALIZATION

	//! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
    //!            in case of restoring from a file with the wrong object type, the file id is printed
    //! @return    class name
    virtual std::string GetTypeId()const
    {
    	return std::string("CallbackHandlerPython");
    }

private:
    PyObject *mCallbackSetParameters;
    PyObject *mCallbackObjective;
    PyObject *mCallbackGradient;
    PyObject *mCallbackHessian;
    
};
} //namespace NuTo
