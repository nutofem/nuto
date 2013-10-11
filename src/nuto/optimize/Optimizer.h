#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include <boost/serialization/vector.hpp>
#include <string>

#include "nuto/base/NuToObject.h"
#include "nuto/optimize/OptimizeException.h"
#include "nuto/optimize/CallbackHandler.h"
#include "nuto/optimize/CallbackHandlerGrid.h"
#include <cfloat>
#include <limits.h>
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
    typedef enum
    {
        MAXFUNCTIONCALLS=1,   //maximum number of function calls is reached
        MAXGRADIENTCALLS,   //maximum number of gradient calls is reached
        MAXHESSIANCALLS,   //maximum number of hessian calls is reached
        MAXITERATIONS,     //maximum number of iterations is reached
        NORMGRADIENT,       //norm of gradient is smaller than a prescribed value
        MINOBJECTIVE,       //objective is smaller than a prescribed value
        DELTAOBJECTIVEBETWEENCYCLES,  //decrease in objective function between two consecutive cycles is smaller than prescribed value
        REACHINGMACHINEPRECISION  //machine precision is reached for the norm of the increment
    } optimization_return_attributes;
    
    Optimizer(unsigned int rNumParameters,unsigned int rNumEqualConstraints,unsigned int rNumInEqualConstraints) : NuToObject()
    {
        mvParameters.Resize(rNumParameters,1);
        mvEqualConstraints.resize(rNumEqualConstraints);
        mvInEqualConstraints.resize(rNumInEqualConstraints);
        mMaxFunctionCalls = INT_MAX;
        mMinObjective = -DBL_MAX;
        isBuild = false;
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
           & BOOST_SERIALIZATION_NVP(objective)
           & BOOST_SERIALIZATION_NVP(mvParameters)
           & BOOST_SERIALIZATION_NVP(mParameters)
           & BOOST_SERIALIZATION_NVP(mvEqualConstraints)
           & BOOST_SERIALIZATION_NVP(mvInEqualConstraints)
           & BOOST_SERIALIZATION_NVP(isBuild)
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
	
    void SetParameters(const NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rParameters)
    {
        mvParameters = rParameters;
        isBuild = false;
    }

    const NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& GetParameters()const
    {
        return mvParameters;    
    }

    void SetParameters(std::vector<double>& rParameters)
    {
    	mParameters=rParameters;
//		assert(mParameters.size()==rParameters.size());
//		for(size_t i=0;i<mParameters.size();++i)
//			mParameters[i]=rParameters[i];
        isBuild = false;
    }

	void GetParameters(std::vector<double>& rParameters)const
	{
		assert(mParameters.size()==rParameters.size());
		for(size_t i=0;i<mParameters.size();++i)
			rParameters[i]=mParameters[i];
	}

    inline int GetNumParameters()
    {
        if(mvParameters.GetNumRows())
        	return mvParameters.GetNumRows();
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
        if (isBuild==true)
            return objective;
        else
        {
            if (mpCallbackHandler!=0)
                return mpCallbackHandler->Objective();
            else
                throw OptimizeException(std::string("[Optimizer::GetObjective]No callback functions defined."));
        }
    }
    
    //! @brief ... Info routine that prints general information about the object (detail according to verbose level)
	virtual void InfoBase()const
	{
		std::cout << "Number of parameters             :" << mvParameters.GetNumRows() << std::endl;
		std::cout << "Number of equality constraints   :" << mvEqualConstraints.size() << std::endl;
		std::cout << "Number of inequality constraints :" << mvInEqualConstraints.size() << std::endl;
		std::cout << "Build                            :" << isBuild << std::endl;
		std::cout << "MaxFunctionCalls                 :" << mMaxFunctionCalls << std::endl;
		std::cout << "MinObjective                     :" << mMinObjective << std::endl;
    }
	
protected:
    CallbackHandler *mpCallbackHandler;
    CallbackHandlerGrid *mpCallbackHandlerGrid;
    double objective;
    FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> mvParameters;
    std::vector<double> mParameters;
    std::vector<double> mvEqualConstraints;
    std::vector<double> mvInEqualConstraints;
    
    bool isBuild;
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
#endif // OPTIMIZER_H
