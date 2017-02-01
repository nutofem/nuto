#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/vector.hpp>
#endif //ENABLE_SERIALIZATION

#include <cfloat>
#include <string>
#include <eigen3/Eigen/Core>

#include "base/NuToObject.h"

#include "optimize/CallbackHandler.h"
#include "optimize/CallbackHandlerGrid.h"
#include "optimize/OptimizeException.h"

namespace NuTo
{


//! @author Joerg F. Unger, ISM
//! @date September 2009
//! @brief ... standard abstract class for all optimizers in NuTo
class Optimizer : public NuToObject
{
#ifdef ENABLE_SERIALIZATION
	friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    enum class eOptimizationReturnAttributes
    {
        MAXFUNCTIONCALLS=1,   //maximum number of function calls is reached
        MAXGRADIENTCALLS,   //maximum number of gradient calls is reached
        MAXHESSIANCALLS,   //maximum number of hessian calls is reached
        MAXITERATIONS,     //maximum number of iterations is reached
        NORMGRADIENT,       //norm of gradient is smaller than a prescribed value
        MINOBJECTIVE,       //mObjective is smaller than a prescribed value
        DELTAOBJECTIVEBETWEENCYCLES,  //decrease in mObjective function between two consecutive cycles is smaller than prescribed value
        REACHINGMACHINEPRECISION  //machine precision is reached for the norm of the increment
    };
    
    Optimizer(unsigned int rNumParameters,unsigned int rNumEqualConstraints,unsigned int rNumInEqualConstraints) : NuToObject()
    {
        mvParameters.resize(rNumParameters,1);
        mvEqualConstraints.resize(rNumEqualConstraints);
        mvInEqualConstraints.resize(rNumInEqualConstraints);
        mMaxFunctionCalls = INT_MAX;
        mMinObjective = -DBL_MAX;
        mIsBuild = false;
	}
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {    
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NuToObject)
           & BOOST_SERIALIZATION_NVP(mpCallbackHandler)
           & BOOST_SERIALIZATION_NVP(mpCallbackHandlerGrid)
           & BOOST_SERIALIZATION_NVP(mObjective)
           & BOOST_SERIALIZATION_NVP(mvParameters)
           & BOOST_SERIALIZATION_NVP(mParameters)
           & BOOST_SERIALIZATION_NVP(mvEqualConstraints)
           & BOOST_SERIALIZATION_NVP(mvInEqualConstraints)
           & BOOST_SERIALIZATION_NVP(mIsBuild)
           & BOOST_SERIALIZATION_NVP(mMinObjective);
    }
#endif // SWIG
#endif // ENABLE_SERIALIZATION

    void SetCallback(NuTo::CallbackHandler* rpCallbackHandler)
	{
	    mpCallbackHandler = rpCallbackHandler;
	}

    void SetCallback(NuTo::CallbackHandlerGrid* rpCallbackHandler)
	{
	    mpCallbackHandlerGrid = rpCallbackHandler;
	}
	
	virtual int Optimize()=0;
	
    void SetParameters(const Eigen::MatrixXd& rParameters)
    {
        mvParameters = rParameters;
        mIsBuild = false;
    }

    const Eigen::MatrixXd& GetParameters()const
    {
        return mvParameters;    
    }

    void SetParameters(std::vector<double>& rParameters)
    {
    	mParameters=rParameters;
        mIsBuild = false;
    }

    std::vector<double>& GetParametersVec()
	{
		return mParameters;
	}

    inline int GetNumParameters()
    {
        if(mvParameters.rows())
        	return mvParameters.rows();
        else
        	return mParameters.size();
    }

    inline void SetMaxFunctionCalls(int rMaxFunctionCalls)
    {
        mMaxFunctionCalls = rMaxFunctionCalls;
    }

    inline void SetMinObjective(double rMinObjective)
    {
        mMinObjective = rMinObjective;
    }

     inline double GetObjective()
    {
        if (mIsBuild)
            return mObjective;
        else
        {
            if (mpCallbackHandler!=0)
                return mpCallbackHandler->Objective();
            else
                throw OptimizeException(__PRETTY_FUNCTION__, "No callback functions defined.");
        }
    }
    
    //! @brief Info routine that prints general information about the object (detail according to verbose level)
    virtual void InfoBase() const
	{
		std::cout << "Number of parameters             :" << mvParameters.rows() << std::endl;
		std::cout << "Number of equality constraints   :" << mvEqualConstraints.size() << std::endl;
		std::cout << "Number of inequality constraints :" << mvInEqualConstraints.size() << std::endl;
		std::cout << "Build                            :" << mIsBuild << std::endl;
		std::cout << "MaxFunctionCalls                 :" << mMaxFunctionCalls << std::endl;
		std::cout << "MinObjective                     :" << mMinObjective << std::endl;
    }
	
protected:
    CallbackHandler *mpCallbackHandler;
    CallbackHandlerGrid *mpCallbackHandlerGrid;
    double mObjective;
    Eigen::MatrixXd mvParameters;
    std::vector<double> mParameters;
    std::vector<double> mvEqualConstraints;
    std::vector<double> mvInEqualConstraints;
    
    bool mIsBuild;
    int mMaxFunctionCalls;
    double mMinObjective;
        
};
} //namespace NuTo
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
#include <boost/serialization/assume_abstract.hpp>
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::Optimize)
#endif // SWIG
#endif // ENABLE_SERIALIZATION

