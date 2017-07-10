#include "mechanics/dofSubMatrixStorage/BlockFullMatrix.h"
#include "base/CallbackInterface.h"
#include "base/Timer.h"
#include "mechanics/structures/StructureOutputBlockMatrix.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/groups/Group.h"
#include "mechanics/structures/StructureBase.h"
#include "mechanics/timeIntegration/ResultGroupNodeForce.h"
#include "mechanics/timeIntegration/TimeIntegrationEnum.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveCalculateStaticData.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveScalar.h"
#include "mechanics/elements/ElementBase.h"
#include "mechanics/structures/StructureBaseEnum.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "math/SparseMatrix.h"
#include "math/SparseMatrixCSRVector2General.h"
#include "mechanics/structures/Assembler.h"

using namespace NuTo;

NewmarkDirect::NewmarkDirect(StructureBase* rStructure)
    : NewmarkBase(rStructure)

{
}


void NewmarkDirect::Info() const
{
    NewmarkBase::Info();
}


void NewmarkDirect::Solve(double rTimeDelta)
{
    Timer timerFull(__FUNCTION__, mStructure->GetShowTime(), mStructure->GetLogger());

    // renumber dofs and build constraint matrix
    mStructure->NodeBuildGlobalDofs(__PRETTY_FUNCTION__);
    const DofStatus& dofStatus = mStructure->GetDofStatus();

    StructureOutputBlockVector delta_dof_dt0(dofStatus, true);

    // [0] = hessian0. , [1] = hessian1 , [2] = hessian2
    std::vector<StructureOutputBlockMatrix> hessian_dt = {StructureOutputBlockMatrix(dofStatus, true),
                                                          StructureOutputBlockMatrix(dofStatus),
                                                          StructureOutputBlockMatrix(dofStatus)};

    // [0] = disp. , [1] = vel. , [2] = acc.
    std::vector<StructureOutputBlockVector> dof_dt = {StructureOutputBlockVector(dofStatus, true),
                                                      StructureOutputBlockVector(dofStatus, true),
                                                      StructureOutputBlockVector(dofStatus, true)};

    // [0] = disp. , [1] = vel. , [2] = acc.
    std::vector<StructureOutputBlockVector> lastConverged_dof_dt = {StructureOutputBlockVector(dofStatus, true),
                                                                    StructureOutputBlockVector(dofStatus, true),
                                                                    StructureOutputBlockVector(dofStatus, true)};

    StructureOutputBlockVector extForce(dofStatus, true);
    StructureOutputBlockVector intForce(dofStatus, true);
    StructureOutputBlockVector residual(dofStatus, true);

    StructureOutputBlockVector prevIntForce(dofStatus, true);
    StructureOutputBlockVector prevExtForce(dofStatus, true);
    StructureOutputBlockVector prevResidual(dofStatus, true);

    BlockFullVector<double> residual_mod(dofStatus);
    BlockFullVector<double> bRHS(dofStatus);
    BlockFullVector<double> deltaBRHS(dofStatus);

    const auto& constraintMatrix = mStructure->GetAssembler().GetConstraintMatrix();

    /*---------------------------------*\
    |         Declare Output Maps       |
    \*---------------------------------*/

    StructureOutputMap evaluate_InternalGradient_Hessian0Hessian1;
    StructureOutputMap evaluate_Hessian0_Hessian1;

    PreIteration(evaluate_InternalGradient_Hessian0Hessian1, evaluate_Hessian0_Hessian1, intForce, hessian_dt,
                 lastConverged_dof_dt, constraintMatrix, residual, residual_mod, dofStatus);

    MainTimeLoop(evaluate_InternalGradient_Hessian0Hessian1, extForce, residual, prevExtForce, deltaBRHS, hessian_dt,
                 lastConverged_dof_dt, residual_mod, constraintMatrix, delta_dof_dt0, dof_dt, prevResidual, rTimeDelta,
                 evaluate_InternalGradient_Hessian0Hessian1, bRHS);
}


std::pair<int, BlockScalar> NewmarkDirect::FindEquilibrium(
        StructureOutputBlockVector& structureResidual,
        StructureOutputMap evalHessian0And1,
        std::vector<StructureOutputBlockMatrix>& hessian_dt, StructureOutputBlockVector& intForce,
        StructureOutputBlockVector& extForce, StructureOutputBlockVector& delta_dof_dt0,
        std::vector<StructureOutputBlockVector>& dof_dt, const BlockSparseMatrix& constraintMatrix, double timeStep)
{
    const auto& dofStatus = mStructure->GetDofStatus();
    BlockFullVector<double> residual(dofStatus);
    structureResidual.ApplyCMatrix(residual, constraintMatrix);
    BlockScalar normResidual = residual.CalculateInfNorm();

    int iteration = 0;
    while (!(normResidual < mToleranceResidual) && iteration < mMaxNumIterations)
    {
        auto inputMap = CreateInputMap();
        mStructure->Evaluate(inputMap, evalHessian0And1);

        if (mCheckCoefficientMatrix)
            mStructure->ElementCheckHessian0(1.e-6, 1.e-8);

        delta_dof_dt0.J = BuildHessianModAndSolveSystem(hessian_dt, residual, timeStep);
        delta_dof_dt0.K = constraintMatrix * delta_dof_dt0.J * (-1.);
        ++mIterationCount;

        double alpha = 1;
        BlockScalar trialNormResidual(mStructure->GetDofStatus());
        StructureOutputBlockVector trial_dof_dt0(dofStatus, true);
        StructureOutputBlockVector trial_dof_dt1(dofStatus, true);
        StructureOutputBlockVector trial_dof_dt2(dofStatus, true);

        // perform a line search
        do
        {
            trial_dof_dt0 = dof_dt[0] + delta_dof_dt0 * alpha;
            if (mStructure->GetNumTimeDerivatives() >= 1)
                trial_dof_dt1 = dof_dt[1] + delta_dof_dt0 * (alpha * mGamma / (timeStep * mBeta));
            if (mStructure->GetNumTimeDerivatives() >= 2)
                trial_dof_dt2 = dof_dt[2] + delta_dof_dt0 * (alpha / (timeStep * timeStep * mBeta));

            MergeDofValues(trial_dof_dt0, trial_dof_dt1, trial_dof_dt2, false);

            auto intForce = EvaluateInternalGradient();

            structureResidual = CalculateResidual(intForce, extForce, hessian_dt[2], trial_dof_dt1, trial_dof_dt2);
            structureResidual.ApplyCMatrix(residual, constraintMatrix);

            trialNormResidual = residual.CalculateInfNorm();

            mStructure->GetLogger() << "[Linesearch a=" << std::to_string(alpha).substr(0, 6)
                                    << "] Trial residual: " << trialNormResidual << "\n";

            alpha *= 0.5;

        } while (mPerformLineSearch && alpha > mMinLineSearchStep && trialNormResidual > (1. - alpha) * normResidual);

        if (alpha > mMinLineSearchStep || !mPerformLineSearch)
        {
            // improvement is achieved, go to next Newton step
            dof_dt[0] = trial_dof_dt0;
            if (mStructure->GetNumTimeDerivatives() >= 1)
                dof_dt[1] = trial_dof_dt1;
            if (mStructure->GetNumTimeDerivatives() >= 2)
                dof_dt[2] = trial_dof_dt2;
            normResidual = trialNormResidual;

            PrintInfoIteration(normResidual, iteration);
            iteration++;
        }
        else
            iteration = mMaxNumIterations;
    }
    return std::make_pair(iteration, normResidual);
}


void NewmarkDirect::MainTimeLoop(
        std::map<NuTo::eStructureOutput, NuTo::StructureOutputBase*>& evaluate_Hessian0_Hessian1,
        StructureOutputBlockVector& extForce, StructureOutputBlockVector& residual,
        StructureOutputBlockVector& prevExtForce, BlockFullVector<double>& deltaBRHS,
        std::vector<StructureOutputBlockMatrix>& hessian_dt,
        std::vector<StructureOutputBlockVector>& lastConverged_dof_dt, BlockFullVector<double>& residual_mod,
        const BlockSparseMatrix& constraintMatrix, StructureOutputBlockVector& delta_dof_dt0,
        std::vector<StructureOutputBlockVector>& dof_dt, StructureOutputBlockVector& prevResidual, double rTimeDelta,
        StructureOutputMap& evaluate_InternalGradient_Hessian0Hessian1, BlockFullVector<double>& bRHS)
{
    const DofStatus& dofStatus = mStructure->GetDofStatus();

    // the minimal time step defined, which is equivalent to six cut-backs
    if (mAutomaticTimeStepping)
    {
        SetMinTimeStep(mMinTimeStep > 0. ? mMinTimeStep : mTimeStep * std::pow(0.5, 6.));
    }
    while (mTimeObject.GetCurrentTime() < rTimeDelta)
    {
        // LEAVE IT OR FIX IT IF YOU CAN:
        // If you dont make a copy of the dof set with dofStatus.GetDofTypes() and use it directly in the for loop
        // everything seems to work fine, but valgrind tells you otherwise. Problem is, that the DofType is changed
        // during the function call inside the loop which leads to reads in already freed blocks of memory
        std::set<Node::eDof> currentDofSet = dofStatus.GetDofTypes();
        for (const auto& dof : currentDofSet)
        {
            mStructure->DofTypeSetIsActive(dof, true);
        }

        if (mTimeObject.GetTimestep() < mMinTimeStep)
            throw MechanicsException(__PRETTY_FUNCTION__,
                                     "time step is smaller than minimum - no convergence is obtained.");

        // calculate Delta_BRhs and Delta_ExtForce
        bRHS = UpdateAndGetAndMergeConstraintRHS(mTimeObject.GetCurrentTime(), lastConverged_dof_dt[0]);
        prevExtForce = CalculateCurrentExternalLoad(mTimeObject.GetCurrentTime());

        mTimeObject.Proceed();

        deltaBRHS = UpdateAndGetConstraintRHS(mTimeObject.GetCurrentTime()) - bRHS;

        IterateForActiveDofValues(evaluate_InternalGradient_Hessian0Hessian1, extForce, residual, prevExtForce,
                                  deltaBRHS, hessian_dt, lastConverged_dof_dt, residual_mod, constraintMatrix,
                                  delta_dof_dt0, dof_dt, prevResidual);

    } // while time steps
}


void NewmarkDirect::CalculateMuDampingMatrix(StructureOutputBlockMatrix& rHessian_dt1,
                                             const StructureOutputBlockMatrix& rHessian_dt2) const
{
    if (mMuDampingMass != 0)
    {
        if (mStructure->GetNumTimeDerivatives() < 2)
            throw MechanicsException(std::string("[") + __PRETTY_FUNCTION__ +
                                     "] MuDampingMass requires a mass matrix (2nd time derivatives).");
        rHessian_dt1.Resize(mStructure->GetDofStatus().GetNumActiveDofsMap(),
                            mStructure->GetDofStatus().GetNumDependentDofsMap());
        rHessian_dt1.SetZero();
        rHessian_dt1.AddScal(rHessian_dt2, mMuDampingMass);
    }
}

StructureOutputBlockVector NewmarkDirect::CalculateDof1(const StructureOutputBlockVector& rDeltaDof_dt0,
                                                        const StructureOutputBlockVector& rDof_dt1,
                                                        const StructureOutputBlockVector& rDof_dt2,
                                                        double rTimeStep) const
{
    return rDeltaDof_dt0 * (mGamma / (rTimeStep * mBeta)) + rDof_dt1 * (1. - mGamma / mBeta) +
           rDof_dt2 * (rTimeStep * (1. - mGamma / (2. * mBeta)));
}


StructureOutputBlockVector NewmarkDirect::CalculateDof2(const StructureOutputBlockVector& rDeltaDof_dt0,
                                                        const StructureOutputBlockVector& rDof_dt1,
                                                        const StructureOutputBlockVector& rDof_dt2,
                                                        double rTimeStep) const
{
    return rDeltaDof_dt0 * (1. / (rTimeStep * rTimeStep * mBeta)) - rDof_dt1 * (1. / (rTimeStep * mBeta)) -
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


void NewmarkDirect::CalculateResidualKForPostprocessing(StructureOutputBlockVector& rResidual,
                                                        const StructureOutputBlockMatrix& rHessian_dt2,
                                                        const StructureOutputBlockVector& rDof_dt1,
                                                        const StructureOutputBlockVector& rDof_dt2) const
{
    if (mStructure->GetDofStatus().HasInteractingConstraints())
        return; // in this case, residual.K is needed for the calculation of residual mod and it is already calculated.

    bool hasNodeForce = false;
    for (const auto& it : mResultMap)
        if (it.second->GetResultType() == eTimeIntegrationResultType::GROUP_NODE_FORCE)
        {
            hasNodeForce = true;
            break; // exit loop
        }
    if (hasNodeForce && mStructure->GetNumTimeDerivatives() >= 2)
    {
        auto dof = rDof_dt1 * mMuDampingMass + rDof_dt2;
        rResidual.K -= rHessian_dt2.KJ * dof.J + rHessian_dt2.KK * dof.K;
    }
    // else:  no need to calculate forces if they are not needed in the post processing
}


void NewmarkDirect::CalculateResidualTrial(StructureOutputBlockVector& rResidual,
                                           const BlockFullVector<double>& rDeltaBRHS,
                                           const std::vector<StructureOutputBlockMatrix>& rHessian_dt,
                                           const StructureOutputBlockVector& rDof_dt1,
                                           const StructureOutputBlockVector& rDof_dt2, double rTimeStep) const
{
    StructureOutputBlockVector deltaDof_dt0(mStructure->GetDofStatus(), true);
    deltaDof_dt0.J.SetZero();
    deltaDof_dt0.K = rDeltaBRHS;

    rResidual -= rHessian_dt[0] * deltaDof_dt0;

    if (mStructure->GetNumTimeDerivatives() >= 1)
    {
        StructureOutputBlockVector delta_dof1 = CalculateDof1(deltaDof_dt0, rDof_dt1, rDof_dt2, rTimeStep) - rDof_dt1;
        rResidual -= rHessian_dt[1] * delta_dof1;
    }

    if (mStructure->GetNumTimeDerivatives() >= 2)
    {
        StructureOutputBlockVector delta_dof2 = CalculateDof2(deltaDof_dt0, rDof_dt1, rDof_dt2, rTimeStep) - rDof_dt2;
        rResidual -= rHessian_dt[2] * delta_dof2;
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


void NewmarkDirect::PrintInfoIteration(const BlockScalar& rNormResidual, int rIteration) const
{
    if (mStructure->GetDofStatus().GetDofTypes().find(Node::eDof::NONLOCALEQSTRAIN) !=
        mStructure->GetDofStatus().GetDofTypes().end())
        return; // Hi Volker. Nothing wrong with the 1 line iteration output in the Solve() method. IMO.
    // And please consider line search...

    Logger& logger = mStructure->GetLogger();

    switch (mVerboseLevel)
    {
    case 0:
    {
        if (rIteration == 0)
        {
            logger << "Iteration:";
        }
        else
        {
            logger << "*";
        }
        return;
    }

    case 1:
    {
        logger << "Iteration: " << rIteration << "\n";
        for (auto dof : mStructure->GetDofStatus().GetActiveDofTypes())
        {
            logger << "Residual " << Node::DofToString(dof) << ": " << rNormResidual[dof] << "\n";
        }
        logger << "--------------------------\n";
        return;
    }

    default:
    {
        return;
    }
    }
}

void NuTo::NewmarkDirect::PreIteration(
        std::map<NuTo::eStructureOutput, NuTo::StructureOutputBase*>& rEvaluate_InternalGradient_Hessian0Hessian1,
        std::map<NuTo::eStructureOutput, NuTo::StructureOutputBase*>& rEvaluate_Hessian0_Hessian1,
        NuTo::StructureOutputBlockVector& rIntForce, std::vector<StructureOutputBlockMatrix>& hessian_dt,
        std::vector<StructureOutputBlockVector>& lastConverged_dof_dt, const BlockSparseMatrix& cmat,
        StructureOutputBlockVector& residual, BlockFullVector<double>& residual_mod, const NuTo::DofStatus& dofStatus)
{
    CalculateStaticAndTimeDependentExternalLoad();

    mToleranceResidual.DefineDefaultValueToIninitializedDofTypes(mToleranceForce);

    if (mMaxTimeStep == 0)
        throw MechanicsException(__PRETTY_FUNCTION__, "max time step is set to zero.");

    if (mStepActiveDofs.empty())
    {
        mStepActiveDofs.push_back(mStructure->DofTypesGetActive());
    }
    else
    {
        for (unsigned int i = 0; i < mStepActiveDofs.size(); ++i)
        {
            if (mStepActiveDofs[i].empty())
            {
                throw MechanicsException(__PRETTY_FUNCTION__,
                                         "Calculation step " + std::to_string(i) + " has no active DOFs.");
            }
        }
    }

    FillOutputMaps(rEvaluate_InternalGradient_Hessian0Hessian1, rEvaluate_Hessian0_Hessian1, rIntForce, hessian_dt,
                   dofStatus);


    ExtractDofValues(lastConverged_dof_dt[0], lastConverged_dof_dt[1], lastConverged_dof_dt[2]);

    UpdateAndGetAndMergeConstraintRHS(mTimeObject.GetCurrentTime(), lastConverged_dof_dt[0]);

    auto inputMap = CreateInputMap();
    mStructure->Evaluate(inputMap, rEvaluate_InternalGradient_Hessian0Hessian1);

    StructureOutputBlockVector initialExtForce =
            CalculateCurrentExternalLoad(mTimeObject.GetCurrentTime()); // put this in element evaluate soon!

    residual = CalculateResidual(rIntForce, initialExtForce, hessian_dt[2], lastConverged_dof_dt[1],
                                 lastConverged_dof_dt[2]);
    residual.ApplyCMatrix(residual_mod, cmat);

    if (mToleranceResidual < residual_mod.CalculateInfNorm())
    {
        mStructure->GetLogger() << residual_mod.CalculateInfNorm();
        throw MechanicsException(__PRETTY_FUNCTION__, "Initial configuration is not in (dynamic) equilibrium.");
    }
    CalculateResidualKForPostprocessing(residual, hessian_dt[2], lastConverged_dof_dt[1], lastConverged_dof_dt[2]);
    PostProcess(residual);
}


void NuTo::NewmarkDirect::FillOutputMaps(
        std::map<NuTo::eStructureOutput, NuTo::StructureOutputBase*>& rEvaluate_InternalGradient_Hessian0Hessian1,
        std::map<NuTo::eStructureOutput, NuTo::StructureOutputBase*>& rEvaluate_Hessian0_Hessian1,
        StructureOutputBlockVector& rIntForce, std::vector<StructureOutputBlockMatrix>& hessian_dt,
        const DofStatus& dofStatus)
{
    rEvaluate_InternalGradient_Hessian0Hessian1[eStructureOutput::INTERNAL_GRADIENT] = &rIntForce;

    rEvaluate_InternalGradient_Hessian0Hessian1[eStructureOutput::HESSIAN0] = &hessian_dt[0];
    rEvaluate_Hessian0_Hessian1[eStructureOutput::HESSIAN0] = &hessian_dt[0];

    if (mStructure->GetNumTimeDerivatives() >= 1 && mMuDampingMass == 0.)
    {
        hessian_dt[1].Resize(dofStatus.GetNumActiveDofsMap(), dofStatus.GetNumDependentDofsMap());
        rEvaluate_InternalGradient_Hessian0Hessian1[eStructureOutput::HESSIAN1] = &hessian_dt[1];
        rEvaluate_Hessian0_Hessian1[eStructureOutput::HESSIAN1] = &hessian_dt[1];
    }

    if (mStructure->GetNumTimeDerivatives() >= 2)
    {
        hessian_dt[2].Resize(dofStatus.GetNumActiveDofsMap(), dofStatus.GetNumDependentDofsMap());
        if (mUseLumpedMass)
        {
            hessian_dt[2] = mStructure->BuildGlobalHessian2Lumped();
        }
        else
        {
            hessian_dt[2] = mStructure->BuildGlobalHessian2();
        }
        CalculateMuDampingMatrix(hessian_dt[1], hessian_dt[2]);
    }
}


ConstitutiveInputMap NuTo::NewmarkDirect::CreateInputMap()
{
    ConstitutiveInputMap inputMap;
    inputMap[Constitutive::eInput::CALCULATE_STATIC_DATA] =
            std::make_unique<ConstitutiveCalculateStaticData>(eCalculateStaticData::EULER_BACKWARD);
    inputMap.emplace(Constitutive::eInput::TIME, std::make_unique<ConstitutiveScalar>());

    double& inputTime = (*inputMap.find(Constitutive::eInput::TIME)->second)[0];
    inputTime = mTimeObject.GetCurrentTime();
    return inputMap;
}


StructureOutputBlockVector NuTo::NewmarkDirect::EvaluateInternalGradient()
{
    auto inputMap = CreateInputMap();
    const auto dofStatus = mStructure->GetDofStatus();
    StructureOutputMap evaluate_InternalGradient;
    StructureOutputBlockVector intForce(dofStatus, true);
    evaluate_InternalGradient[eStructureOutput::INTERNAL_GRADIENT] = &intForce;
    mStructure->Evaluate(inputMap, evaluate_InternalGradient);
    return intForce;
}


void NewmarkDirect::IterateForActiveDofValues(
        std::map<NuTo::eStructureOutput, NuTo::StructureOutputBase*>& evaluate_Hessian0_Hessian1,
        StructureOutputBlockVector& extForce, StructureOutputBlockVector& residual,
        StructureOutputBlockVector& prevExtForce, BlockFullVector<double>& deltaBRHS,
        std::vector<StructureOutputBlockMatrix>& hessian_dt,
        std::vector<StructureOutputBlockVector>& lastConverged_dof_dt, BlockFullVector<double>& residual_mod,
        const BlockSparseMatrix& constraintMatrix, StructureOutputBlockVector& delta_dof_dt0,
        std::vector<StructureOutputBlockVector>& dof_dt,
        StructureOutputBlockVector& prevResidual)
{
    // at the moment needed to do the postprocessing after the last step and not after every step of a staggered
    // solution.
    unsigned int staggeredStepNumber = 0;
    for (const auto& activeDofs : mStepActiveDofs)
    {
        ++staggeredStepNumber;
        mStructure->DofTypeSetIsActive(activeDofs);

        PrintInfoStagger();

        auto inputMap = CreateInputMap();
        mStructure->Evaluate(inputMap, evaluate_Hessian0_Hessian1);

        extForce = CalculateCurrentExternalLoad(mTimeObject.GetCurrentTime());

        /*------------------------------------------------*\
        |         Calculate Residual for trail state       |
        \*------------------------------------------------*/
        residual = extForce - prevExtForce;
        CalculateResidualTrial(residual, deltaBRHS, hessian_dt, lastConverged_dof_dt[1], lastConverged_dof_dt[2],
                               mTimeObject.GetTimestep());
        residual.ApplyCMatrix(residual_mod, constraintMatrix);

        mStructure->GetLogger() << "\n"
                                << "Initial trial residual:               " << residual_mod.CalculateInfNorm() << "\n";

        delta_dof_dt0.J = BuildHessianModAndSolveSystem(hessian_dt, residual_mod, mTimeObject.GetTimestep());
        delta_dof_dt0.K = deltaBRHS - constraintMatrix * delta_dof_dt0.J;
        ++mIterationCount;

        // calculate trial state
        dof_dt[0] = lastConverged_dof_dt[0] + delta_dof_dt0;
        if (mStructure->GetNumTimeDerivatives() >= 1)
            dof_dt[1] = CalculateDof1(delta_dof_dt0, lastConverged_dof_dt[1], lastConverged_dof_dt[2],
                                      mTimeObject.GetTimestep());
        if (mStructure->GetNumTimeDerivatives() >= 2)
            dof_dt[2] = CalculateDof2(delta_dof_dt0, lastConverged_dof_dt[1], lastConverged_dof_dt[2],
                                      mTimeObject.GetTimestep());


        MergeDofValues(dof_dt[0], dof_dt[1], dof_dt[2], false);

        auto intForce = EvaluateInternalGradient();

        residual = CalculateResidual(intForce, extForce, hessian_dt[2], dof_dt[1], dof_dt[2]);

        std::pair<int, BlockScalar> result =
                FindEquilibrium(residual, evaluate_Hessian0_Hessian1, hessian_dt,
                                intForce, extForce, delta_dof_dt0, dof_dt, constraintMatrix, mTimeObject.GetTimestep());

        auto iterations = result.first;
        auto residualNorm = result.second;

        if (residualNorm < mToleranceResidual)
        {
            // converged solution
            if (mVerboseLevel > 2)
                mStructure->GetLogger() << "\n *** UpdateStaticData *** from NewMarkDirect \n";

            // Update static data
            mStructure->ElementTotalUpdateStaticData();

            // store converged step
            lastConverged_dof_dt[0] = dof_dt[0];
            if (mStructure->GetNumTimeDerivatives() >= 1)
                lastConverged_dof_dt[1] = dof_dt[1];
            if (mStructure->GetNumTimeDerivatives() >= 1)
                lastConverged_dof_dt[2] = dof_dt[2];

            prevResidual = residual;

            MergeDofValues(dof_dt[0], dof_dt[1], dof_dt[2], true);

            mTime += mTimeObject.GetTimestep();

            mStructure->GetLogger() << "Convergence after " << iterations << " iterations at time " << mTime
                                    << " (timestep " << mTimeObject.GetTimestep() << ").\n";
            mStructure->GetLogger() << "Residual: \t" << residualNorm << "\n";
            // perform Postprocessing
            if (staggeredStepNumber >= mStepActiveDofs.size())
            {
                CalculateResidualKForPostprocessing(prevResidual, hessian_dt[2], dof_dt[1], dof_dt[2]);
                PostProcess(prevResidual);
            }

            // eventually increase next time step
            if (mAutomaticTimeStepping && iterations < 0.25 * mMaxNumIterations)
            {
                mTimeObject.ScaleTimestep(1.5);
            }

            if (mCallback && mCallback->Exit(*mStructure))
                return;
        }
        else
        {
            mStructure->GetLogger() << "No convergence with timestep " << mTimeObject.GetTimestep() << " at time "
                                    << mTime << "\n";
            // no convergence
            if (mAutomaticTimeStepping)
            {
                // no convergence, reduce the time step and start from scratch
                mTimeObject.RestorePreviosTime();
                mTimeObject.ScaleTimestep(0.5);
            }
            else
            {
                throw MechanicsException(__PRETTY_FUNCTION__, "No convergence with the current maximum number of "
                                                              "iterations, either use automatic time stepping, "
                                                              "reduce the time step or the minimal line search cut "
                                                              "back factor.");
            }
        } // if tolerance
    } // active dof loop
}


BlockFullVector<double>
NewmarkDirect::BuildHessianModAndSolveSystem(std::vector<StructureOutputBlockMatrix>& rHessian_dt,
                                             const BlockFullVector<double>& rResidualMod, double rTimeStep) const
{
    Timer timer(__FUNCTION__, GetShowTime(), mStructure->GetLogger());

    // since rHessian0 will change in the next iteration, the rHessian0 will be the hessian for the solver
    if (mStructure->GetNumTimeDerivatives() >= 1)
        rHessian_dt[0].AddScal(rHessian_dt[1], mGamma / (mBeta * rTimeStep));

    if (mStructure->GetNumTimeDerivatives() >= 2)
        rHessian_dt[0].AddScal(rHessian_dt[2], 1. / (mBeta * rTimeStep * rTimeStep));

    rHessian_dt[0].ApplyCMatrix(mStructure->GetAssembler().GetConstraintMatrix());

    auto result = mSolver->Solve(rHessian_dt[0].JJ, rResidualMod);
    return result;
}
