#pragma once

#include "mechanics/timeIntegration/NewmarkBase.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "mechanics/structures/StructureBaseEnum.h"

namespace NuTo
{

class StructureOutputBlockMatrix;
class StructureOutputBlockVector;
enum class eStructureOutput;
namespace Constitutive
{
    enum class eInput;

}
template <typename IOEnum> class ConstitutiveIOMap;
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

    //! @brief ... Info routine that prints general information about the object (detail according to verbose level)
    void Info() const override;

protected:
    using StructureOutputMap = std::map<eStructureOutput, StructureOutputBase*>;

    StructureOutputBlockVector CalculateDof1(const StructureOutputBlockVector& rDeltaDof_dt0,
                                             const StructureOutputBlockVector& rDof_dt1,
                                             const StructureOutputBlockVector& rDof_dt2, double rTimeStep) const;

    StructureOutputBlockVector CalculateDof2(const StructureOutputBlockVector& rDeltaDof_dt0,
                                             const StructureOutputBlockVector& rDof_dt1,
                                             const StructureOutputBlockVector& rDof_dt2, double rTimeStep) const;

    //! @brief ... builds the modified hessian matrix (including cmat) and solves the system
    //! @param rHessian_dt0 ... hessian_dt0 matrix
    //! @param rHessian_dt1 ... hessian_dt1 matrix
    //! @param rHessian_dt2 ... hessian_dt2 matrix
    //! @param rResidualMod ... modified residual (including the cmatrix)
    //! @param rTimeStep ... current time step
    //! @return ... deltaDof_dt0.J
    //! @remark ... If hessian_dt0 is constant, its values are preserved (hessianMod = temporary matrix). Otherwise,
    //! hessian0 will be used (hessianMod = hessian0)
    BlockFullVector<double> BuildHessianModAndSolveSystem(StructureOutputBlockMatrix& rHessian_dt0,
                                                          const StructureOutputBlockMatrix& rHessian_dt1,
                                                          const StructureOutputBlockMatrix& rHessian_dt2,
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

    void CalculateMuDampingMatrix(StructureOutputBlockMatrix& rHessian_dt1,
                                  const StructureOutputBlockMatrix& rHessian_dt2) const;

    void CalculateResidualTrial(StructureOutputBlockVector& rResidual, const BlockFullVector<double>& rDeltaBRHS,
                                const StructureOutputBlockMatrix& rHessian_dt0,
                                const StructureOutputBlockMatrix& rHessian_dt1,
                                const StructureOutputBlockMatrix& rHessian_dt2,
                                const StructureOutputBlockVector& rDof_dt1, const StructureOutputBlockVector& rDof_dt2,
                                double rTimeStep) const;


    //! @brief Prints Info about the current calculation stage
    void PrintInfoStagger() const;

    //! @brief Prints Info about the current iteration
    void PrintInfoIteration(const BlockScalar& rNormResidual, int rIteration) const;


    // temporary stuff (no doxy-doc)

    void PreIteration(std::map<NuTo::eStructureOutput, NuTo::StructureOutputBase *> &rEvaluate_InternalGradient,
                      std::map<NuTo::eStructureOutput, NuTo::StructureOutputBase *> &rEvaluate_InternalGradient_Hessian0Hessian1,
                      std::map<NuTo::eStructureOutput, NuTo::StructureOutputBase *> &rEvaluate_Hessian0_Hessian1,
                      NuTo::ConstitutiveInputMap &rInputMap,
                      NuTo::StructureOutputBlockVector &rIntForce,
                      std::vector<NuTo::StructureOutputBlockMatrix> &hessianVec, NuTo::StructureOutputBlockVector &lastConverged_dof_dt0, NuTo::StructureOutputBlockVector &lastConverged_dof_dt1, NuTo::StructureOutputBlockVector &lastConverged_dof_dt2, const NuTo::BlockSparseMatrix &cmat, NuTo::StructureOutputBlockVector &residual, BlockFullVector<double> &residual_mod, double curTime,
                      const NuTo::DofStatus &dofStatus);

    void FillOutputMaps(std::map<NuTo::eStructureOutput, NuTo::StructureOutputBase *> &rEvaluate_InternalGradient,
                        std::map<NuTo::eStructureOutput, NuTo::StructureOutputBase *> &rEvaluate_InternalGradient_Hessian0Hessian1,
                        std::map<NuTo::eStructureOutput, NuTo::StructureOutputBase *> &rEvaluate_Hessian0_Hessian1,
                        NuTo::StructureOutputBlockVector &rIntForce,
                        std::vector<NuTo::StructureOutputBlockMatrix>& hessianVec,
                        const NuTo::DofStatus &dofStatus);

    void FillInputMap(NuTo::ConstitutiveInputMap &rInputMap);

    void IterateForActiveDofValues(NuTo::ConstitutiveInputMap& inputMap,
                                   std::map<NuTo::eStructureOutput,
                                   NuTo::StructureOutputBase *> &evaluate_Hessian0_Hessian1,
                                   NuTo::StructureOutputBlockVector &extForce,
                                   double &curTime,
                                   NuTo::StructureOutputBlockVector &residual,
                                   NuTo::StructureOutputBlockVector &prevExtForce,
                                   BlockFullVector<double> &deltaBRHS,
                                   std::vector<StructureOutputBlockMatrix>& hessianVec,
                                   NuTo::StructureOutputBlockVector &lastConverged_dof_dt0,
                                   NuTo::StructureOutputBlockVector &lastConverged_dof_dt1,
                                   NuTo::StructureOutputBlockVector &lastConverged_dof_dt2,
                                   double &timeStep,
                                   BlockFullVector<double> &residual_mod,
                                   const NuTo::BlockSparseMatrix &constraintMatrix,
                                   NuTo::StructureOutputBlockVector &delta_dof_dt0,
                                   NuTo::StructureOutputBlockVector &dof_dt0,
                                   NuTo::StructureOutputBlockVector &dof_dt1,
                                   NuTo::StructureOutputBlockVector &dof_dt2,
                                   StructureOutputMap &evaluate_InternalGradient,
                                   NuTo::StructureOutputBlockVector &intForce,
                                   NuTo::StructureOutputBlockVector &prevResidual,
                                   double& inputTime);


    std::pair<int, BlockScalar> FindEquilibrium(StructureOutputBlockVector& structureResidual, const ConstitutiveInputMap& inputMap,
                                StructureOutputMap evalHessian0And1, StructureOutputMap evalGradient,
                                std::vector<NuTo::StructureOutputBlockMatrix>& hessianVec,
                                StructureOutputBlockVector& intForce,
                                StructureOutputBlockVector& extForce,
                                StructureOutputBlockVector& delta_dof_dt0,
                                StructureOutputBlockVector& dof_dt0, StructureOutputBlockVector& dof_dt1,
                                StructureOutputBlockVector& dof_dt2, const BlockSparseMatrix& constraintMatrix,
                                double timeStep);

    void MainTimeLoop(NuTo::ConstitutiveInputMap &inputMap,
                      std::map<NuTo::eStructureOutput,
                      NuTo::StructureOutputBase *> &evaluate_Hessian0_Hessian1,
                      NuTo::StructureOutputBlockVector &extForce,
                      double &curTime,
                      NuTo::StructureOutputBlockVector &residual,
                      NuTo::StructureOutputBlockVector &prevExtForce,
                      BlockFullVector<double> &deltaBRHS,
                      std::vector<NuTo::StructureOutputBlockMatrix>& hessianVec,
                      NuTo::StructureOutputBlockVector &lastConverged_dof_dt0,
                      NuTo::StructureOutputBlockVector &lastConverged_dof_dt1,
                      NuTo::StructureOutputBlockVector &lastConverged_dof_dt2,
                      double &timeStep,
                      BlockFullVector<double> &residual_mod,
                      const NuTo::BlockSparseMatrix &constraintMatrix,
                      NuTo::StructureOutputBlockVector &delta_dof_dt0,
                      NuTo::StructureOutputBlockVector &dof_dt0,
                      NuTo::StructureOutputBlockVector &dof_dt1,
                      NuTo::StructureOutputBlockVector &dof_dt2,
                      StructureOutputMap &evaluate_InternalGradient,
                      NuTo::StructureOutputBlockVector &intForce,
                      NuTo::StructureOutputBlockVector &prevResidual,
                      double rTimeDelta,
                      StructureOutputMap &evaluate_InternalGradient_Hessian0Hessian1,
                      BlockFullVector<double> &bRHS);

protected:
    double mMinLineSearchStep = 0.01;
    bool mPerformLineSearch = true;

    int mVerboseLevel = 1; //!< controls the output verbosity (0 = silent)
};
} // namespace NuTo
