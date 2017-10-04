#include "base/CallbackInterface.h"
#include "base/Timer.h"
#include "mechanics/structures/StructureOutputBlockMatrix.h"
#include "mechanics/structures/StructureBase.h"
#include "mechanics/timeIntegration/postProcessing/ResultGroupNodeForce.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveCalculateStaticData.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveScalar.h"
#include "mechanics/structures/StructureBaseEnum.h"
#include "mechanics/timeIntegration/postProcessing/PostProcessor.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "math/SparseMatrix.h"
#include "math/SparseMatrixCSRVector2General.h"
#include "mechanics/structures/Assembler.h"
#include "mechanics/nodes/NodeEnum.h"

using namespace NuTo;

NewmarkDirect::NewmarkDirect(StructureBase* rStructure)
    : TimeIntegrationBase(rStructure)
    , mHessian2(rStructure->GetDofStatus())
    , mDampingMatrix(rStructure->GetDofStatus())
{
}


void NewmarkDirect::Solve(double timeFinal)
{
    mTimeControl.SetTimeFinal(timeFinal);
    Timer timerFull(__FUNCTION__, mStructure->GetShowTime(), mStructure->GetLogger());

    mStructure->NodeBuildGlobalDofs(__PRETTY_FUNCTION__);
    auto dofValues = InitialState();

    if (mAutomaticTimeStepping && mTimeControl.GetMinTimeStep() <= 0.)
    {
        // the minimal time step is equivalent to six cut-backs
        mTimeControl.SetMinTimeStep(mTimeControl.GetTimeStep() * std::pow(0.5, 6.));
    }


    while (not mTimeControl.Finished())
    {
        mTimeControl.Proceed();

        // LEAVE IT OR FIX IT IF YOU CAN:
        // If you dont make a copy of the dof set with dofStatus.GetDofTypes() and use it directly in the for loop
        // everything seems to work fine, but valgrind tells you otherwise. Problem is, that the DofType is changed
        // during the function call inside the loop which leads to reads in already freed blocks of memory
        std::set<Node::eDof> currentDofSet = mStructure->GetDofStatus().GetDofTypes();
        for (const auto& dof : currentDofSet)
        {
            mStructure->DofTypeSetIsActive(dof, true);
        }

        const auto bRHS = UpdateAndGetAndMergeConstraintRHS(mTimeControl.GetPreviousTime(), dofValues[0]);
        const auto prevExtForce = CalculateCurrentExternalLoad(mTimeControl.GetPreviousTime());
        const auto deltaBRHS = UpdateAndGetConstraintRHS(mTimeControl.GetCurrentTime()) - bRHS;

        IterateForActiveDofValues(prevExtForce, deltaBRHS, dofValues);
    }
}


std::pair<int, BlockScalar> NewmarkDirect::FindEquilibrium(StructureOutputBlockVector& structureResidual,
                                                           const StructureOutputBlockVector& extForce,
                                                           StructureOutputBlockVector& delta_dof_dt0,
                                                           std::array<StructureOutputBlockVector, 3>& dof_dt)
{
    const auto& dofStatus = mStructure->GetDofStatus();
    const auto& constraintMatrix = mStructure->GetAssembler().GetConstraintMatrix();
    auto residual = Assembler::ApplyCMatrix(structureResidual, constraintMatrix);
    BlockScalar residualNorm = residual.CalculateInfNorm();

    int iteration = 0;
    while (not(residualNorm < mToleranceResidual) and iteration < mMaxNumIterations)
    {
        auto hessians = EvaluateHessians();

        delta_dof_dt0.J = BuildHessianModAndSolveSystem(hessians, residual, mTimeControl.GetTimeStep());

        delta_dof_dt0.K = constraintMatrix * delta_dof_dt0.J * (-1.);
        ++mIterationCount;

        double alpha = 1;
        BlockScalar trialNormResidual(dofStatus);
        StructureOutputBlockVector trial_dof_dt0(dofStatus, true);
        StructureOutputBlockVector trial_dof_dt1(dofStatus, true);
        StructureOutputBlockVector trial_dof_dt2(dofStatus, true);

        // perform a line search
        do
        {
            trial_dof_dt0 = dof_dt[0] + delta_dof_dt0 * alpha;
            if (mStructure->GetNumTimeDerivatives() >= 1)
                trial_dof_dt1 = dof_dt[1] + delta_dof_dt0 * (alpha * mGamma / (mTimeControl.GetTimeStep() * mBeta));
            if (mStructure->GetNumTimeDerivatives() >= 2)
                trial_dof_dt2 =
                        dof_dt[2] +
                        delta_dof_dt0 * (alpha / (mTimeControl.GetTimeStep() * mTimeControl.GetTimeStep() * mBeta));

            MergeDofValues(trial_dof_dt0, trial_dof_dt1, trial_dof_dt2, false);

            const auto intForce = EvaluateInternalGradient();

            structureResidual = CalculateResidual(intForce, extForce, hessians[2], trial_dof_dt1, trial_dof_dt2);
            residual = Assembler::ApplyCMatrix(structureResidual, constraintMatrix);

            trialNormResidual = residual.CalculateInfNorm();

            mStructure->GetLogger() << "[Linesearch a=" << std::to_string(alpha).substr(0, 6)
                                    << "] Trial residual: " << trialNormResidual << "\n";

            alpha *= 0.5;

        } while (mPerformLineSearch && alpha > mMinLineSearchStep && trialNormResidual > (1. - alpha) * residualNorm);

        if (alpha > mMinLineSearchStep || !mPerformLineSearch)
        {
            // improvement is achieved, go to next Newton step
            dof_dt[0] = trial_dof_dt0;
            if (mStructure->GetNumTimeDerivatives() >= 1)
                dof_dt[1] = trial_dof_dt1;
            if (mStructure->GetNumTimeDerivatives() >= 2)
                dof_dt[2] = trial_dof_dt2;
            residualNorm = trialNormResidual;

            PrintInfoIteration(residualNorm, iteration);
            iteration++;
        }
        else
            iteration = mMaxNumIterations;
    }
    return std::make_pair(iteration, residualNorm);
}


StructureOutputBlockMatrix NewmarkDirect::CalculateMuDampingMatrix(const StructureOutputBlockMatrix& hessian2) const
{
    if (mUseMuDamping and mStructure->GetNumTimeDerivatives() < 2)
        throw Exception(__PRETTY_FUNCTION__, "MuDampingMass requires a mass matrix (2nd time derivatives).");

    const auto& dofStatus = mStructure->GetDofStatus();
    StructureOutputBlockMatrix hessian1(dofStatus);
    hessian1.Resize(dofStatus.GetNumActiveDofsMap(), dofStatus.GetNumDependentDofsMap());
    hessian1.SetZero();
    hessian1.AddScal(hessian2, mMuDampingMass);
    return hessian1;
}

StructureOutputBlockVector NewmarkDirect::CalculateDof1(const StructureOutputBlockVector& rDeltaDof_dt0,
                                                        const StructureOutputBlockVector& rDof_dt1,
                                                        const StructureOutputBlockVector& rDof_dt2) const
{
    const auto step = mTimeControl.GetTimeStep();
    return rDeltaDof_dt0 * (mGamma / (step * mBeta)) + rDof_dt1 * (1. - mGamma / mBeta) +
           rDof_dt2 * (step * (1. - mGamma / (2. * mBeta)));
}


StructureOutputBlockVector NewmarkDirect::CalculateDof2(const StructureOutputBlockVector& rDeltaDof_dt0,
                                                        const StructureOutputBlockVector& rDof_dt1,
                                                        const StructureOutputBlockVector& rDof_dt2) const
{
    const auto step = mTimeControl.GetTimeStep();
    return rDeltaDof_dt0 * (1. / (step * step * mBeta)) - rDof_dt1 * (1. / (step * mBeta)) -
           rDof_dt2 * ((0.5 - mBeta) / mBeta);
}


StructureOutputBlockVector NewmarkDirect::CalculateResidual(const StructureOutputBlockVector& rIntForce,
                                                            const StructureOutputBlockVector& rExtForce,
                                                            const StructureOutputBlockMatrix& rHessian_dt2,
                                                            const StructureOutputBlockVector& rDof_dt1,
                                                            const StructureOutputBlockVector& rDof_dt2) const
{
    StructureOutputBlockVector residual = rExtForce - rIntForce;

    // The residual for numTimeDerivatives = 1 is included in the internal forces.
    // If there is muDamping, there must be a rHessian2.
    if (mStructure->GetNumTimeDerivatives() >= 2)
    {
        residual -= rHessian_dt2 * (rDof_dt1 * mMuDampingMass + rDof_dt2);
    }
    return residual;
}


void NewmarkDirect::CalculateResidualTrial(StructureOutputBlockVector& rResidual,
                                           const BlockFullVector<double>& rDeltaBRHS,
                                           const std::array<StructureOutputBlockMatrix, 3>& hessians,
                                           const StructureOutputBlockVector& rDof_dt1,
                                           const StructureOutputBlockVector& rDof_dt2) const
{
    StructureOutputBlockVector deltaDof_dt0(mStructure->GetDofStatus(), true);
    deltaDof_dt0.J.SetZero();
    deltaDof_dt0.K = rDeltaBRHS;

    rResidual -= hessians[0] * deltaDof_dt0;

    if (mStructure->GetNumTimeDerivatives() >= 1)
    {
        StructureOutputBlockVector delta_dof1 = CalculateDof1(deltaDof_dt0, rDof_dt1, rDof_dt2) - rDof_dt1;
        rResidual -= hessians[1] * delta_dof1;
    }

    if (mStructure->GetNumTimeDerivatives() >= 2)
    {
        StructureOutputBlockVector delta_dof2 = CalculateDof2(deltaDof_dt0, rDof_dt1, rDof_dt2) - rDof_dt2;
        rResidual -= hessians[2] * delta_dof2;
    }
}


void NewmarkDirect::PrintInfoStagger() const
{
    if (mStepActiveDofs.size() == 1) // equals unstaggered solution, no info needed
        return;

    Logger& logger = mStructure->GetLogger();
    logger << "\n"
           << "Activated Dofs: | ";
    for (auto dof : mStructure->GetDofStatus().GetActiveDofTypes())
    {
        logger << Node::DofToString(dof) << " | ";
    }
    logger << "\n";
}


void NewmarkDirect::PrintInfoIteration(const BlockScalar& residualNorm, int iteration) const
{
    Logger& logger = mStructure->GetLogger();
    switch (mVerboseLevel)
    {
    case 0:
    {
        if (iteration == 0)
        {
            logger << "Iteration:";
        }
        else
        {
            logger << "*";
        }
        break;
    }
    case 1:
    {
        logger << "Iteration: " << iteration << "\n";
        for (const auto dof : mStructure->GetDofStatus().GetActiveDofTypes())
        {
            logger << "Residual " << Node::DofToString(dof) << ": " << residualNorm[dof] << "\n";
        }
        logger << "--------------------------\n";
        break;
    }
    default:
        break;
    }
}


std::array<StructureOutputBlockVector, 3> NuTo::NewmarkDirect::InitialState()
{
    CalculateStaticAndTimeDependentExternalLoad();

    mToleranceResidual.DefineDefaultValueToIninitializedDofTypes(mToleranceForce);

    if (mStepActiveDofs.empty())
        mStepActiveDofs.push_back(mStructure->DofTypesGetActive());
    else
    {
        for (const auto& activeDofs : mStepActiveDofs)
        {
            if (activeDofs.empty())
                throw Exception(__PRETTY_FUNCTION__, "Calculation step has no active DOFs.");
        }
    }

    if (mStructure->GetNumTimeDerivatives() >= 2)
    {
        if (mUseLumpedMass)
            mHessian2 = mStructure->BuildGlobalHessian2Lumped();
        else
            mHessian2 = mStructure->BuildGlobalHessian2();

        mDampingMatrix = CalculateMuDampingMatrix(mHessian2);
    }

    auto dofValues = ExtractDofValues();

    UpdateAndGetAndMergeConstraintRHS(mTimeControl.GetCurrentTime(), dofValues[0]);

    const auto gradientAndHessians = EvaluateGradientAndHessians();
    const auto intForce = gradientAndHessians.first;
    const auto hessians = gradientAndHessians.second;

    const auto initialExtForce = CalculateCurrentExternalLoad(mTimeControl.GetCurrentTime());

    auto residual = CalculateResidual(intForce, initialExtForce, hessians[2], dofValues[1], dofValues[2]);
    const auto& constraintMatrix = mStructure->GetAssembler().GetConstraintMatrix();
    const auto residual_mod = Assembler::ApplyCMatrix(residual, constraintMatrix);

//    if (mToleranceResidual < residual_mod.CalculateInfNorm())
//    {
//        mStructure->GetLogger() << residual_mod.CalculateInfNorm();
//        throw Exception(__PRETTY_FUNCTION__, "Initial configuration is not in (dynamic) equilibrium.");
//    }

    mPostProcessor->PostProcess(residual);
    return dofValues;
}


ConstitutiveInputMap NuTo::NewmarkDirect::CreateInputMap()
{
    ConstitutiveInputMap inputMap;
    inputMap[Constitutive::eInput::CALCULATE_STATIC_DATA] =
            std::make_unique<ConstitutiveCalculateStaticData>(eCalculateStaticData::EULER_BACKWARD);
    inputMap.emplace(Constitutive::eInput::TIME, std::make_unique<ConstitutiveScalar>());

    double& inputTime = (*inputMap.find(Constitutive::eInput::TIME)->second)[0];
    inputTime = mTimeControl.GetCurrentTime();
    return inputMap;
}


StructureOutputBlockVector NuTo::NewmarkDirect::EvaluateInternalGradient()
{
    const auto& dofStatus = mStructure->GetDofStatus();
    StructureOutputMap evalInternalGradient;
    StructureOutputBlockVector intForce(dofStatus, true);
    evalInternalGradient[eStructureOutput::INTERNAL_GRADIENT] = &intForce;

    auto inputMap = CreateInputMap();
    mStructure->Evaluate(inputMap, evalInternalGradient);
    return intForce;
}


std::array<StructureOutputBlockMatrix, 3> NuTo::NewmarkDirect::EvaluateHessians()
{
    const auto& dofStatus = mStructure->GetDofStatus();
    std::array<StructureOutputBlockMatrix, 3> hessians = {StructureOutputBlockMatrix(dofStatus, true),
                                                          StructureOutputBlockMatrix(dofStatus), mHessian2};
    StructureOutputMap evalHessians;
    evalHessians[eStructureOutput::HESSIAN0] = &hessians[0];
    if (mStructure->GetNumTimeDerivatives() >= 1)
    {
        if (mUseMuDamping)
            hessians[1] = mDampingMatrix;
        else
        {
            hessians[1].Resize(dofStatus.GetNumActiveDofsMap(), dofStatus.GetNumDependentDofsMap());
            evalHessians[eStructureOutput::HESSIAN1] = &hessians[1];
        }
    }

    auto inputMap = CreateInputMap();
    mStructure->Evaluate(inputMap, evalHessians);
    if (mCheckCoefficientMatrix)
        mStructure->ElementCheckHessian0(1.e-6, 1.e-8);

    return hessians;
}


std::pair<StructureOutputBlockVector, std::array<StructureOutputBlockMatrix, 3>>
NuTo::NewmarkDirect::EvaluateGradientAndHessians()
{
    const auto& dofStatus = mStructure->GetDofStatus();
    StructureOutputBlockVector intForce(dofStatus, true);
    std::array<StructureOutputBlockMatrix, 3> hessians = {StructureOutputBlockMatrix(dofStatus, true),
                                                          StructureOutputBlockMatrix(dofStatus), mHessian2};
    StructureOutputMap evalGradAndHessians;
    evalGradAndHessians[eStructureOutput::HESSIAN0] = &hessians[0];
    evalGradAndHessians[eStructureOutput::INTERNAL_GRADIENT] = &intForce;
    if (mStructure->GetNumTimeDerivatives() >= 1)
    {
        if (mUseMuDamping)
            hessians[1] = mDampingMatrix;
        else
        {
            hessians[1].Resize(dofStatus.GetNumActiveDofsMap(), dofStatus.GetNumDependentDofsMap());
            evalGradAndHessians[eStructureOutput::HESSIAN1] = &hessians[1];
        }
    }

    auto inputMap = CreateInputMap();
    mStructure->Evaluate(inputMap, evalGradAndHessians);
    return std::make_pair(intForce, hessians);
}


void NewmarkDirect::IterateForActiveDofValues(const StructureOutputBlockVector& prevExtForce,
                                              const BlockFullVector<double>& deltaBRHS,
                                              std::array<StructureOutputBlockVector, 3>& lastConverged_dof_dt)
{
    // at the moment needed to do the postprocessing after the last step and not after every step of a staggered
    // solution.
    unsigned int staggeredStepNumber = 0;
    const auto& dofStatus = mStructure->GetDofStatus();
    // [0] = disp. , [1] = vel. , [2] = acc.
    std::array<StructureOutputBlockVector, 3> dof_dt = {StructureOutputBlockVector(dofStatus, true),
                                                        StructureOutputBlockVector(dofStatus, true),
                                                        StructureOutputBlockVector(dofStatus, true)};
    StructureOutputBlockVector delta_dof_dt0(dofStatus, true);
    const auto& constraintMatrix = mStructure->GetAssembler().GetConstraintMatrix();

    // Needed for mTimeControl.Proceed() function
    int timeStepMaxIterations = 0;
    bool converged = true;

    //auto final_lastConverged_dof_dt = lastConverged_dof_dt;

    for (const auto& activeDofs : mStepActiveDofs)
    {
        ++staggeredStepNumber;
        mStructure->DofTypeSetIsActive(activeDofs);

        PrintInfoStagger();

        auto hessians = EvaluateHessians();

        const auto extForce = CalculateCurrentExternalLoad(mTimeControl.GetCurrentTime());

        auto residual = extForce - prevExtForce;
        CalculateResidualTrial(residual, deltaBRHS, hessians, lastConverged_dof_dt[1], lastConverged_dof_dt[2]);
        const auto residual_mod = Assembler::ApplyCMatrix(residual, constraintMatrix);

        mStructure->GetLogger() << "\n"
                                << "Initial trial residual:               " << residual_mod.CalculateInfNorm() << "\n";

        delta_dof_dt0.J = BuildHessianModAndSolveSystem(hessians, residual_mod, mTimeControl.GetTimeStep());

        delta_dof_dt0.K = deltaBRHS - constraintMatrix * delta_dof_dt0.J;
        ++mIterationCount;

        // calculate trial state
        dof_dt[0] = lastConverged_dof_dt[0] + delta_dof_dt0;
        if (mStructure->GetNumTimeDerivatives() >= 1)
            dof_dt[1] = CalculateDof1(delta_dof_dt0, lastConverged_dof_dt[1], lastConverged_dof_dt[2]);
        if (mStructure->GetNumTimeDerivatives() >= 2)
            dof_dt[2] = CalculateDof2(delta_dof_dt0, lastConverged_dof_dt[1], lastConverged_dof_dt[2]);

        MergeDofValues(dof_dt[0], dof_dt[1], dof_dt[2], false);

        auto intForce = EvaluateInternalGradient();
        residual = CalculateResidual(intForce, extForce, hessians[2], dof_dt[1], dof_dt[2]);

        const auto result = FindEquilibrium(residual, extForce, delta_dof_dt0, dof_dt);
        const auto iterations = result.first;
        const auto residualNorm = result.second;

        if (iterations > timeStepMaxIterations)
            timeStepMaxIterations = iterations;

        converged = residualNorm < mToleranceResidual;
        if (converged)
        {
            mStructure->ElementTotalUpdateStaticData();

            auto prevResidual = residual;

            MergeDofValues(dof_dt[0], dof_dt[1], dof_dt[2], true); // only merges currently active dof types

            mStructure->GetLogger() << "Convergence after " << iterations << " iterations at time "
                                    << mTimeControl.GetCurrentTime() << " (timestep " << mTimeControl.GetTimeStep()
                                    << ").\n";
            mStructure->GetLogger() << "Residual: \t" << residualNorm << "\n";

            if (staggeredStepNumber >= mStepActiveDofs.size())
                mPostProcessor->PostProcess(prevResidual);

            if (mCallback && mCallback->Exit(*mStructure))
                return;
        }
        else
        {
            mStructure->GetLogger() << "No convergence with timestep " << mTimeControl.GetTimeStep() << " at time "
                                    << mTimeControl.GetCurrentTime() << "\n";
            break;
        } // if tolerance
    } // active dof loop

    // Continue with next timestep or reduce timestep and restart iteration
    mTimeControl.AdjustTimestep(timeStepMaxIterations, mMaxNumIterations, converged);

    // store converged dofs after each staggered step is done
    if (converged)
        lastConverged_dof_dt = ExtractDofValues();
}


BlockFullVector<double>
        NewmarkDirect::BuildHessianModAndSolveSystem(std::array<StructureOutputBlockMatrix, 3>& rHessians,
                                                     const BlockFullVector<double>& rResidualMod,
                                                     double rTimeStep) const
{
    Timer timer(__FUNCTION__, GetShowTime(), mStructure->GetLogger());

    // since rHessian0 will change in the next iteration, the rHessian0 will be the hessian for the solver
    if (mStructure->GetNumTimeDerivatives() >= 1)
        rHessians[0].AddScal(rHessians[1], mGamma / (mBeta * rTimeStep));

    if (mStructure->GetNumTimeDerivatives() >= 2)
        rHessians[0].AddScal(rHessians[2], 1. / (mBeta * rTimeStep * rTimeStep));

    rHessians[0].ApplyCMatrix(mStructure->GetAssembler().GetConstraintMatrix());

    return mSolver->Solve(rHessians[0].JJ, rResidualMod);
}


void NewmarkDirect::MergeDofValues(const StructureOutputBlockVector& rDof_dt0,
                                   const StructureOutputBlockVector& rDof_dt1,
                                   const StructureOutputBlockVector& rDof_dt2, bool rMergeAll)
{
    mStructure->NodeMergeDofValues(0, rDof_dt0.J, rDof_dt0.K);

    if (mStructure->GetNumTimeDerivatives() >= 1)
    {
        if (rMergeAll or mMuDampingMass == 0)
        {
            mStructure->NodeMergeDofValues(1, rDof_dt1.J, rDof_dt1.K);
        }
    }

    if (mStructure->GetNumTimeDerivatives() >= 2)
    {
        if (rMergeAll)
        {
            mStructure->NodeMergeDofValues(2, rDof_dt2.J, rDof_dt2.K);
        }
    }
    mStructure->ElementTotalUpdateTmpStaticData();
}
