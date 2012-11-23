// $Id$
#ifndef CALLBACKHANDLERGRID_H
#define CALLBACKHANDLERGRID_H

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/base/NuToObject.h"
#include "nuto/optimize/OptimizeException.h"
#include <iostream>
#include <vector>

namespace NuTo
{
//! @author Andrea Keszler, ISM
//! @date July 2010
//! @brief ... abstract class to handle callback routines
class CallbackHandlerGrid : public virtual NuToObject
{
#ifdef ENABLE_SERIALIZATION
	friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
    CallbackHandlerGrid(): NuToObject()
    {
    	mCallbackSetParameters = 0;
    	mCallbackGradient = 0;
    	mCallbackHessian = 0;
    }


#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize CallbackHandler" << std::endl;
#endif
    	ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NuToObject)
    	& BOOST_SERIALIZATION_NVP(mCallbackSetParameters)
    	& BOOST_SERIALIZATION_NVP(mCallbackGradient)
    	& BOOST_SERIALIZATION_NVP(mCallbackHessian);
    	mCallbackSetParameters = 0;
    	mCallbackGradient = 0;
    	mCallbackHessian = 0;
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize CallbackHandler" << std::endl;
#endif

	}
#endif // ENABLE_SERIALIZATION

    virtual std::vector<double>&  GetParameters()
    {
		throw OptimizeException("[CallbackHandler::SetParameters] SetParameters function not implemented in CallbackHandler object.");
    }

    virtual std::vector<double>&  GetResidual()
    {
		throw OptimizeException("[CallbackHandler::SetParameters] GetResidual function not implemented in CallbackHandler object.");
    }
    virtual void SetParameters(std::vector<double>& rParameters)
    {
		throw OptimizeException("[CallbackHandler::SetParameters] SetParameters function not implemented in CallbackHandler object.");
    }

    virtual void SetResidual(std::vector<double>& rResidual)
    {
		throw OptimizeException("[CallbackHandler::SetParameters] SetResidual function not implemented in CallbackHandler object.");
    }

    virtual void Gradient (std::vector<double>& rValue,std::vector<double>& rGradient)const
	{
		throw OptimizeException("[CallbackHandler::Gradient] Gradient function not implemented in CallbackHandler object.");
	}

	virtual void Hessian (std::vector<double>& rHessian)const
	{
		throw OptimizeException("[CallbackHandler::Gradient] Gradient function not implemented in CallbackHandler object.");
	}

	void Info()const
	{
		std::cout << "CallbackHandlerGrid" << std::endl;
	}

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

    virtual std::string GetTypeId()const
    {
    	return std::string("CallbackHandlerGrid");
    }

private:
    const std::vector<double> *mCallbackSetParameters;
    std::vector<double> *mCallbackGradient;
    std::vector<double> *mCallbackHessian;
};
} //namespace NuTo
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
#include <boost/serialization/assume_abstract.hpp>
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::CallbackHandlerGrid)
#endif // SWIG
#endif // ENABLE_SERIALIZATION

#endif // CALLBACKHANDLERPYTHON_H
