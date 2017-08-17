#pragma once

#include "mechanics/timeIntegration/TimeIntegrationBase.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "mechanics/structures/StructureBaseEnum.h"
#include "mechanics/structures/StructureOutputBlockMatrix.h"


namespace NuTo
{

class StructureOutputBlockMatrix;
class StructureOutputBlockVector;
enum class eStructureOutput;
namespace Constitutive
{
enum class eInput;
}
template <typename IOEnum>
class ConstitutiveIOMap;
typedef ConstitutiveIOMap<Constitutive::eInput> ConstitutiveInputMap;

//! @brief Newmark time integration scheme
class NewmarkDirect : public TimeIntegrationBase
{

public:
    //! @brief constructor
    NewmarkDirect(StructureBase* rStructure);

    void SetPerformLineSearch(bool rPerformLineSearch)
    {
        mPerformLineSearch = rPerformLineSearch;
    }

    int GetVerboseLevel() const
    {
        return mVerboseLevel;
    }

    void SetVerboseLevel(int rVerboseLevel)
    {
        mVerboseLevel = rVerboseLevel;
    }


    //! @brief returns true, if the method is only conditionally stable (for unconditional stable, this is false)
    bool HasCriticalTimeStep() const override
    {
        return false;
    }

    //! @brief calculate the critical time step for explicit routines
    //! for implicit routines, this will simply return zero (cmp HasCriticalTimeStep())
    double CalculateCriticalTimeStep() const override
    {
        return 0;
    }

    virtual void SetTimeAndTimeStep(double& curTime, double& timeStep, double rTimeDelta)
    {
        // calculate time and time step close to the end time of the integration interval
        if (rTimeDelta - curTime < 0.5 * timeStep)
        {
            timeStep += rTimeDelta - curTime;
            curTime = rTimeDelta;
        }
    }


    //! @brief perform the time integration
    //! @param timeFinal ... simulation continues until time hits this value
    virtual void Solve(double timeFinal) override;

    void SetDampingCoefficientMass(double rMuDampingMass)
    {
        mUseMuDamping = true;
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
    using StructureOutputMap = std::map<eStructureOutput, StructureOutputBase*>;

    StructureOutputBlockVector CalculateDof1(const StructureOutputBlockVector& rDeltaDof_dt0,
                                             const StructureOutputBlockVector& rDof_dt1,
                                             const StructureOutputBlockVector& rDof_dt2) const;

    StructureOutputBlockVector CalculateDof2(const StructureOutputBlockVector& rDeltaDof_dt0,
                                             const StructureOutputBlockVector& rDof_dt1,
                                             const StructureOutputBlockVector& rDof_dt2) const;

    //! @brief ... builds the modified hessian matrix (including cmat) and solves the system
    BlockFullVector<double> BuildHessianModAndSolveSystem(std::array<NuTo::StructureOutputBlockMatrix, 3>& rHessians,
                                                          const BlockFullVector<double>& rResidualMod,
                                                          double rTimeStep) const;

    StructureOutputBlockVector CalculateResidual(const StructureOutputBlockVector& rIntForce,
                                                 const StructureOutputBlockVector& rExtForce,
                                                 const StructureOutputBlockMatrix& rHessian2,
                                                 const StructureOutputBlockVector& rDof_dt1,
                                                 const StructureOutputBlockVector& rDof_dt2) const;

    StructureOutputBlockMatrix CalculateMuDampingMatrix(const StructureOutputBlockMatrix& hessian2) const;

    void CalculateResidualTrial(StructureOutputBlockVector& rResidual, const BlockFullVector<double>& rDeltaBRHS,
                                const std::array<NuTo::StructureOutputBlockMatrix, 3>& hessians,
                                const StructureOutputBlockVector& rDof_dt1,
                                const StructureOutputBlockVector& rDof_dt2) const;


    //! @brief Prints Info about the current calculation stage
    void PrintInfoStagger() const;

    //! @brief Prints Info about the current iteration
    void PrintInfoIteration(const BlockScalar& rNormResidual, int rIteration) const;

private:
    std::array<StructureOutputBlockVector, 3> InitialState();

    void IterateForActiveDofValues(const StructureOutputBlockVector& prevExtForce,
                                   const BlockFullVector<double>& deltaBRHS,
                                   std::array<NuTo::StructureOutputBlockVector, 3>& lastConverged_dof_dt);


    std::pair<int, BlockScalar> FindEquilibrium(StructureOutputBlockVector& structureResidual,
                                                const StructureOutputBlockVector& extForce,
                                                StructureOutputBlockVector& delta_dof_dt0,
                                                std::array<StructureOutputBlockVector, 3>& dof_dt);

    ConstitutiveInputMap CreateInputMap();
    StructureOutputBlockVector EvaluateInternalGradient();
    std::array<StructureOutputBlockMatrix, 3> EvaluateHessians();
    std::pair<StructureOutputBlockVector, std::array<StructureOutputBlockMatrix, 3>> EvaluateGradientAndHessians();


protected:
    double mMinLineSearchStep = 0.01;
    bool mPerformLineSearch = true;

    int mVerboseLevel = 1; //!< controls the output verbosity (0 = silent)


    StructureOutputBlockMatrix mHessian2;

    bool mUseMuDamping = false;
    double mMuDampingMass = 0; //!< damping coefficient for the mass (F^d = -mMuDampingMass*M*v)
    StructureOutputBlockMatrix mDampingMatrix;

    // NewtonRaphson parameters
    double mToleranceForce = 1.e-6;
    int mMaxNumIterations = 20;

    // Newmark parameters
    double mBeta = 0.25;
    double mGamma = 0.5;

    bool mUseLumpedMass = false;
};
} // namespace NuTo
