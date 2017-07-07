#pragma once

#include "mechanics/timeIntegration/TimeIntegrationBase.h"


namespace NuTo
{
//! @author Jörg F. Unger, NU
//! @date February 2012
//! @brief ... standard class for implicit timeintegration (Newmark, but you can use it for statics as well with setting
//! the flag isDynamic to false)
class NewmarkBase : public TimeIntegrationBase
{

public:
    //! @brief constructor
    NewmarkBase(StructureBase* rStructure);

    virtual void Solve(double rTimeDelta) override = 0;

    void SetDampingCoefficientMass(double rMuDampingMass)
    {
        mMuDampingMass = rMuDampingMass;
    }

    void SetToleranceForce(double rToleranceForce)
    {
        mToleranceForce = rToleranceForce;
    }

    void SetMaxNumIterations(int rMaxNumIterations)
    {
        mMaxNumIterations = rMaxNumIterations;
    }

    //! @brief merges the dof values depending on the numTimeDerivatives and rMergeAll
    //! @param rDof_dt0 ... 0th time derivative
    //! @param rDof_dt1 ... 1st time derivative
    //! @param rDof_dt2 ... 2nd time derivative
    //! @param rMergeAll ... false: merges dof_dt1 only when mMuMassDamping = 0, ignores dof_dt2
    void MergeDofValues(const StructureOutputBlockVector& rDof_dt0, const StructureOutputBlockVector& rDof_dt1,
                        const StructureOutputBlockVector& rDof_dt2, bool rMergeAll);


protected:
    double mMuDampingMass = 0; //!< damping coefficient for the mass (F^d = -mMuDampingMass*M*v)

    // NewtonRaphson parameters
    double mToleranceForce = 1.e-6;
    int mMaxNumIterations = 20;

    // Newmark parameters
    double mBeta = 0.25;
    double mGamma = 0.5;

    double mInternalEnergy = 0.;
    double mExternalEnergy = 0.;
    double mKineticEnergy = 0.;
    double mDampedEnergy = 0.;

    bool mUseLumpedMass = false;
};
} // namespace NuTo
