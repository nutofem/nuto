// $Id$

#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/timeIntegration/TimeIntegrationBase.h"

namespace NuTo
{
//! @author Jörg F. Unger, NU
//! @date February 2012
//! @brief ... standard class for implicit timeintegration (Newmark, but you can use it for statics as well with setting the flag isDynamic to false)
class NewmarkBase : public TimeIntegrationBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:

    //! @brief constructor
    NewmarkBase(StructureBase* rStructure);

    virtual NuTo::eError Solve(double rTimeDelta)=0;

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

    bool GetUseLumpedMass()const
    {
    	return mUseLumpedMass;
    }

    void SetUseLumpedMass(bool rUseLumpedMass)
    {
    	mUseLumpedMass = rUseLumpedMass;
    }

    //! @brief merges the dof values depending on the numTimeDerivatives and rMergeAll
    //! @param rDof_dt0 ... 0th time derivative
    //! @param rDof_dt1 ... 1st time derivative
    //! @param rDof_dt2 ... 2nd time derivative
    //! @param rMergeAll ... false: merges dof_dt1 only when mMuMassDamping = 0, ignores dof_dt2
    void MergeDofValues(const StructureOutputBlockVector& rDof_dt0, const StructureOutputBlockVector& rDof_dt1, const StructureOutputBlockVector& rDof_dt2, bool rMergeAll);


#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif// SWIG
#endif // ENABLE_SERIALIZATION

    //! @brief ... Info routine that prints general information about the object (detail according to verbose level)
    void Info()const;


protected:

#ifdef ENABLE_SERIALIZATION
    NewmarkBase(){};
#endif  // ENABLE_SERIALIZATION


    //damping coefficient for the mass (F^d = -mMuDampingMass*M*v)
	double mMuDampingMass;
    //NewtonRaphson parameters
    double mToleranceForce;
	int mMaxNumIterations;
	//Newmark parameters
	double mBeta;
	double mGamma;
	double mInternalEnergy;
	double mExternalEnergy;
	double mKineticEnergy;
	double mDampedEnergy;

	bool mUseLumpedMass;
};
} //namespace NuTo
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_CLASS_EXPORT_KEY(NuTo::NewmarkBase)
#endif // SWIG
#endif // ENABLE_SERIALIZATION



