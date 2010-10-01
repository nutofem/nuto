#ifndef CALLBACKHANDLERGRID_H
#define CALLBACKHANDLERGRID_H

#include "nuto/optimize/CallbackHandler.h"
#include "nuto/optimize/OptimizeException.h"

namespace NuTo
{
//! @author Andrea Keszler, ISM
//! @date July 2010
//! @brief ... abstract class to handle callback routines
class CallbackHandlerGrid : public CallbackHandler
{
#ifdef ENABLE_SERIALIZATION
	friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
    CallbackHandlerGrid();

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
    	ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(CallbackHandler)
    	& BOOST_SERIALIZATION_NVP(mCallbackSetParameters)
    	& BOOST_SERIALIZATION_NVP(mCallbackObjective)
    	& BOOST_SERIALIZATION_NVP(mCallbackGradient)
    	& BOOST_SERIALIZATION_NVP(mCallbackHessian);
;
    	mCallbackSetParameters = 0;
    	mCallbackObjective = 0;
    	mCallbackGradient = 0;
    	mCallbackHessian = 0;

	}
#endif // ENABLE_SERIALIZATION

	void SetParameters(const NuTo::FullMatrix<double>& rParameters);
	void SetParameters(NuTo::FullMatrix<double>& rParameters);

	double Objective()const;

	void Gradient (NuTo::FullMatrix<double>& rGradient)const;

    void Hessian(NuTo::FullMatrix<double>& rHessian)const;

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

    virtual std::string GetTypeId()const
    {
    	return std::string("CallbackHandlerGrid");
    }

private:
    FullMatrix<double> *mCallbackSetParameters;
    double *mCallbackObjective;
    FullMatrix<double> *mCallbackGradient;
    FullMatrix<double> *mCallbackHessian;
};
} //namespace NuTo
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
#include <boost/serialization/assume_abstract.hpp>
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::CallbackHandlerGrid)
#endif // SWIG
#endif // ENABLE_SERIALIZATION

#endif // CALLBACKHANDLERPYTHON_H
