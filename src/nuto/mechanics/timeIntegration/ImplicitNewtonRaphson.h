// $Id$

#ifndef IMPLICIT_NEWTON_RAPHSON_H
#define IMPLICIT_NEWTON_RAPHSON_H

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/timeIntegration/TimeIntegrationBase.h"

namespace NuTo
{
//! @author Jörg F. Unger, NU
//! @date February 2012
//! @brief ... standard class for implicit timeintegration (static Newton Raphson or Newmark for dynamics)
class ImplicitNewtonRaphson : public TimeIntegrationBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:

    //! @brief constructor
    ImplicitNewtonRaphson();

    void SetDampingCoefficientMass(double rMuDampingMass)
    {
    	mMuDampingMass = rMuDampingMass;
    }

    double GetDampingCoefficientMass()const
    {
    	return mMuDampingMass;
    }

    void SetToleranceForce(double rToleranceForce)
    {
    	mToleranceForce = rToleranceForce;
    }

    double GetToleranceForce()const
    {
    	return mToleranceForce;
    }

    void SetMaxNumIterations(int rMaxNumIterations)
    {
    	mMaxNumIterations = rMaxNumIterations;
    }

    int GetMaxNumIterations()const
    {
    	return mMaxNumIterations;
    }

    void SetMinLineSearchStep(double rMinLineSearchStep)
    {
    	mMinLineSearchStep = rMinLineSearchStep;
    }

    double GetMinLineSearchStep()const
    {
    	return mMinLineSearchStep;
    }

    void SetNewmarkBeta(double rBeta)
    {
    	mBeta = rBeta;
    }

    double GetNewmarkBeta()const
    {
    	return mBeta;
    }

    void SetNewmarkGamma(double rGamma)
    {
    	mGamma = rGamma;
    }

    double GetNewmarkGamma()const
    {
    	return mGamma;
    }

    virtual bool IsDynamic()const=0;


#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif// SWIG
#endif // ENABLE_SERIALIZATION

    //! @brief perform the time integration
    //! @param rStructure ... structure
    //! @param rTimeDelta ... length of the simulation
    NuTo::Error::eError Solve(StructureBase& rStructure, double rTimeDelta);

    //! @brief ... Info routine that prints general information about the object (detail according to verbose level)
    void Info()const;

protected:
    //damping coefficient for the mass (F^d = -mMuDampingMass*M*v)
	double mMuDampingMass;
    //NewtonRaphson parameters
	double mToleranceForce;
	int mMaxNumIterations;
	double mMinLineSearchStep;
	//Newmark parameters
	double mBeta;
	double mGamma;
	double mInternalEnergy;
	double mExternalEnergy;
	double mKineticEnergy;
	double mDampedEnergy;
};
} //namespace NuTo
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_CLASS_EXPORT_KEY(NuTo::ImplicitNewtonRaphson)
#endif // SWIG
#endif // ENABLE_SERIALIZATION


#endif // IMPLICIT_NEWTON_RAPHSON_H
