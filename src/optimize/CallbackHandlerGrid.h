#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include <boost/dynamic_bitset.hpp>
#include "optimize/OptimizeException.h"
#include <iostream>
#include <vector>

namespace NuTo
{
//! @author Andrea Keszler, ISM
//! @date July 2010
//! @brief ... abstract class to handle callback routines
class CallbackHandlerGrid
{
#ifdef ENABLE_SERIALIZATION
	friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
    CallbackHandlerGrid()
    {
    	mCallbackSetParameters = 0;
    	mCallbackGradient = 0;
    	mCallbackHessian = 0;
    }

    virtual ~CallbackHandlerGrid() = default;


#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize CallbackHandlerGrid" << std::endl;
#endif
    	ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NuToObject)
    	& BOOST_SERIALIZATION_NVP(mCallbackSetParameters)
    	& BOOST_SERIALIZATION_NVP(mCallbackGradient)
    	& BOOST_SERIALIZATION_NVP(mCallbackHessian);
    	mCallbackSetParameters = 0;
    	mCallbackGradient = 0;
    	mCallbackHessian = 0;
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize CallbackHandlerGrid" << std::endl;
#endif

	}
#endif // ENABLE_SERIALIZATION

    virtual std::vector<double>&  GetParameters()
    {
		throw OptimizeException("[CallbackHandlerGrid::SetParameters] SetParameters function not implemented in CallbackHandlerGrid object.");
    }

    virtual std::vector<double>&  GetRightHandSide()
    {
		throw OptimizeException("[CallbackHandlerGrid::GetRightHandSide] GetRightHandSide function not implemented in CallbackHandlerGrid object.");
    }
    virtual void SetParameters(std::vector<double>& rParameters)
    {
		throw OptimizeException("[CallbackHandlerGrid::SetParameters] SetParameters function not implemented in CallbackHandlerGrid object.");
    }

    virtual void SetRightHandSide(std::vector<double>& rRightHandSide)
    {
		throw OptimizeException("[CallbackHandlerGrid::SetRightHandSide] SetRightHandSide function not implemented in CallbackHandlerGrid object.");
    }

    virtual void Gradient (std::vector<double>& rValue,std::vector<double>& rGradient)
	{
		throw OptimizeException("[CallbackHandlerGrid::Gradient] Gradient function not implemented in CallbackHandlerGrid object.");
	}

 	virtual void Hessian (std::vector<double>& rHessian)
	{
		throw OptimizeException("[CallbackHandlerGrid::Gradient] Gradient function not implemented in CallbackHandlerGrid object.");
	}

	//! @brief get DisplacementConstaints
	//! @return dynamic_bitset of constraints
	virtual const boost::dynamic_bitset<> GetDisplacementConstaints()
	{
		throw OptimizeException("[CallbackHandlerGrid::GetDisplacementConstaints] GetDisplacementConstaints function not implemented in CallbackHandlerGrid object.");
	}

	//! @brief correct solution for hanging nodes
	//! @param displacement solution
	virtual void HangingNodesCorrection(std::vector<double>& u)
	{
		throw OptimizeException("[CallbackHandlerGrid::HangingNodesCorrection] HangingNodesCorrection function not implemented in CallbackHandlerGrid object.");
	}


	virtual void SetMisesWielandt (bool rMisesWielandt)
	{
		throw OptimizeException("[CallbackHandlerGrid::SetMisesWielandt] SetMisesWielandt function not implemented in CallbackHandlerGrid object.");
	}
	virtual double GetWeightingFactor()
	{
		throw OptimizeException("[CallbackHandlerGrid::GetWeightingFactor] GetWeightingFactor function not implemented in CallbackHandlerGrid object.");
	}
	virtual void SetWeightingFactor(double rWeight)
	{
		throw OptimizeException("[CallbackHandlerGrid::GetWeightingFactor] GetWeightingFactor function not implemented in CallbackHandlerGrid object.");

	}
	void Info() const
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

