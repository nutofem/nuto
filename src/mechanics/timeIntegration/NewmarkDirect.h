#pragma once

#include "mechanics/timeIntegration/NewmarkBase.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "mechanics/structures/StructureBaseEnum.h"
#include "mechanics/timeIntegration/Time.h"

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
class NewmarkDirect : public NewmarkBase
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
    //! @param rTimeDelta ... length of the simulation
    virtual void Solve(double rTimeDelta) override;

protected:
    using StructureOutputMap = std::map<eStructureOutput, StructureOutputBase*>;

    StructureOutputBlockVector CalculateDof1(const StructureOutputBlockVector& rDeltaDof_dt0,
                                             const StructureOutputBlockVector& rDof_dt1,
                                             const StructureOutputBlockVector& rDof_dt2, double rTimeStep) const;

    StructureOutputBlockVector CalculateDof2(const StructureOutputBlockVector& rDeltaDof_dt0,
                                             const StructureOutputBlockVector& rDof_dt1,
                                             const StructureOutputBlockVector& rDof_dt2, double rTimeStep) const;

    //! @brief ... builds the modified hessian matrix (including cmat) and solves the system
    BlockFullVector<double> BuildHessianModAndSolveSystem(std::array<NuTo::StructureOutputBlockMatrix, 3>& rHessians,
                                                          const BlockFullVector<double>& rResidualMod,
                                                          double rTimeStep) const;

    StructureOutputBlockVector CalculateResidual(const StructureOutputBlockVector& rIntForce,
                                                 const StructureOutputBlockVector& rExtForce,
                                                 const StructureOutputBlockMatrix& rHessian2,
                                                 const StructureOutputBlockVector& rDof_dt1,
                                                 const StructureOutputBlockVector& rDof_dt2) const;

    //! @brief Calculates (if needed) the residual.K part for the post-processing. Since it is not needed for the actual
    //! time integration
    //! its calculation is skipped if Cmat has only zero entries.
    void CalculateResidualKForPostprocessing(StructureOutputBlockVector& rResidual,
                                             const StructureOutputBlockMatrix& rHessian_dt2,
                                             const StructureOutputBlockVector& rDof_dt1,
                                             const StructureOutputBlockVector& rDof_dt2) const;

    StructureOutputBlockMatrix CalculateMuDampingMatrix(const StructureOutputBlockMatrix& hessian2) const;

    void CalculateResidualTrial(StructureOutputBlockVector& rResidual, const BlockFullVector<double>& rDeltaBRHS,
                                const std::array<NuTo::StructureOutputBlockMatrix, 3>& hessians,
                                const StructureOutputBlockVector& rDof_dt1, const StructureOutputBlockVector& rDof_dt2,
                                double rTimeStep) const;


    //! @brief Prints Info about the current calculation stage
    void PrintInfoStagger() const;

    //! @brief Prints Info about the current iteration
    void PrintInfoIteration(const BlockScalar& rNormResidual, int rIteration) const;

public:
    // workaround until dynamic timestepping is also enabled by time class
    virtual void SetTimeStep(double rTimeStep) override
    {
        mTimeStep = rTimeStep;
        mTimeObject.SetEquidistantTimestepping(rTimeStep);
    }

    //! @brief sets the maximum time step for the time integration procedure
    virtual void SetMaxTimeStep(double rMaxTimeStep) override
    {
        mMaxTimeStep = rMaxTimeStep;
        mTimeObject.SetMaxTimestep(rMaxTimeStep);
    }

    //! @brief sets the minimum time step for the time integration procedure
    virtual void SetMinTimeStep(double rMinTimeStep) override
    {
        mMinTimeStep = rMinTimeStep;
        mTimeObject.SetMinTimestep(rMinTimeStep);
    }

protected:
    std::array<StructureOutputBlockVector, 3> InitialState();

    void IterateForActiveDofValues(const StructureOutputBlockVector& prevExtForce, const BlockFullVector<double>& deltaBRHS,
                                   std::array<NuTo::StructureOutputBlockVector, 3>& lastConverged_dof_dt);


    std::pair<int, BlockScalar> FindEquilibrium(StructureOutputBlockVector& structureResidual,
                                                const StructureOutputBlockVector& extForce,
                                                StructureOutputBlockVector& delta_dof_dt0,
                                                std::array<StructureOutputBlockVector, 3>& dof_dt,
                                                const BlockSparseMatrix& constraintMatrix, double timeStep);

    ConstitutiveInputMap CreateInputMap();
    StructureOutputBlockVector EvaluateInternalGradient();
    std::array<StructureOutputBlockMatrix, 3> EvaluateHessians();
    std::pair<StructureOutputBlockVector, std::array<StructureOutputBlockMatrix, 3>> EvaluateGradientAndHessians();

protected:
    double mMinLineSearchStep = 0.01;
    bool mPerformLineSearch = true;

    int mVerboseLevel = 1; //!< controls the output verbosity (0 = silent)

    Time mTimeObject;
    StructureOutputBlockMatrix mHessian2;
};
} // namespace NuTo
