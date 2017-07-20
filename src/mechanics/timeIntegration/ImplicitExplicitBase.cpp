/*
 * ImplicitExplicitBase.cpp
 *
 *  Created on: 17 Mar 2016
 *      Author: ttitsche
 */

#include "base/CallbackInterface.h"
#include "math/SparseMatrixCSRGeneral_Def.h"
#include "mechanics/timeIntegration/ImplicitExplicitBase.h"
#include "base/Timer.h"
#include "mechanics/structures/StructureBase.h"
#include "mechanics/structures/StructureBaseEnum.h"
#include "mechanics/structures/StructureOutputBlockMatrix.h"
#include "mechanics/structures/StructureOutputDummy.h"
#include "math/SparseDirectSolverMUMPS.h"

#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveTimeStep.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveCalculateStaticData.h"
#include "mechanics/structures/Assembler.h"
#include "mechanics/timeIntegration/postProcessing/PostProcessor.h"

NuTo::ImplicitExplicitBase::ImplicitExplicitBase(StructureBase* rStructure)
    : TimeIntegrationBase(rStructure)
{
}

NuTo::ImplicitExplicitBase::~ImplicitExplicitBase()
{
}

void NuTo::ImplicitExplicitBase::Solve(double rTimeDelta)
{
    NuTo::Timer timerFull(__FUNCTION__, mStructure->GetShowTime(), mStructure->GetLogger());
    NuTo::Timer timerDebug("Init", true, mStructure->GetLogger());

    mStructure->NodeBuildGlobalDofs(__PRETTY_FUNCTION__);

    CalculateStaticAndTimeDependentExternalLoad();

    const DofStatus& dofStatus = mStructure->GetDofStatus();

    if (mStepActiveDofs.empty())
        throw Exception(__PRETTY_FUNCTION__, "Define a set of active dofs for each calculation step. ");

    // deactivate all dof types
    mStructure->DofTypeDeactivateAll();

    std::map<Node::eDof, SparseDirectSolverMUMPS> preFactorizedHessians;
    FactorizeConstantHessians(preFactorizedHessians);


    /*---------------------------------*\
    |        Allocate Variables         |
    \*---------------------------------*/

    StructureOutputBlockMatrix hessian0(dofStatus, true);

    StructureOutputBlockVector delta_dof_dt0(dofStatus, true);

    StructureOutputBlockVector dof_dt0(dofStatus, true); // e.g. disp
    StructureOutputBlockVector lastConverged_dof_dt0(dofStatus, true); // e.g. disp

    StructureOutputBlockVector extForce(dofStatus, true);
    StructureOutputBlockVector intForce(dofStatus, true);
    StructureOutputBlockVector residual(dofStatus, true);


    // for constraints
    // ---------------

    BlockFullVector<double> residual_mod(dofStatus);

    /*---------------------------------*\
    |    Declare and fill Output Maps   |
    \*---------------------------------*/

    // Declare output maps
    //        std::map<NuTo::eStructureOutput, NuTo::StructureOutputBase*> evalUpdateStaticData;
    StructureOutputDummy dummy;
    //        evalUpdateStaticData                [eStructureOutput::UPDATE_STATIC_DATA] = &dummy;

    std::map<NuTo::eStructureOutput, NuTo::StructureOutputBase*> evalInternalGradient;
    evalInternalGradient[eStructureOutput::INTERNAL_GRADIENT] = &intForce;
    evalInternalGradient[eStructureOutput::UPDATE_STATIC_DATA] = &dummy;

    std::map<NuTo::eStructureOutput, NuTo::StructureOutputBase*> evalInternalGradientAndHessian0;
    evalInternalGradientAndHessian0[eStructureOutput::INTERNAL_GRADIENT] = &intForce;
    evalInternalGradientAndHessian0[eStructureOutput::HESSIAN0] = &hessian0;
    evalInternalGradientAndHessian0[eStructureOutput::UPDATE_STATIC_DATA] = &dummy;

    /*---------------------------------*\
    |    Declare and fill Input map     |
    \*---------------------------------*/


    ConstitutiveInputMap input;
    input[Constitutive::eInput::TIME_STEP] = std::make_unique<ConstitutiveTimeStep>(2);
    auto& timeStep = *static_cast<ConstitutiveTimeStep*>(input[Constitutive::eInput::TIME_STEP].get());
    input[Constitutive::eInput::CALCULATE_STATIC_DATA] =
            std::make_unique<ConstitutiveCalculateStaticData>(eCalculateStaticData::USE_PREVIOUS, 1);
    auto& calculateStaticData =
            *static_cast<ConstitutiveCalculateStaticData*>(input[Constitutive::eInput::CALCULATE_STATIC_DATA].get());

    timeStep.SetCurrentTimeStep(mTimeStep);
    timeStep.SetCurrentTimeStep(mTimeStep);
    timeStep.SetCurrentTimeStep(mTimeStep);
    timeStep.SetCurrentTimeStep(mTimeStep);
    timeStep.SetCurrentTimeStep(mTimeStep);

    mPostProcessor->PostProcess(residual);

    while (mTime < rTimeDelta)
    {
        timerDebug.Reset("\033[1;31m Iteration " + std::to_string(mIterationCount) + " at current time: " +
                         std::to_string(mTime) + "\033[0m");

        mStructure->DofTypeActivateAll();


        // calculate Delta_BRhs and Delta_ExtForce
        auto bRHS = UpdateAndGetConstraintRHS(mTime);
        auto prevExtForce = CalculateCurrentExternalLoad(mTime);

        mTime += mTimeStep;

        auto deltaBRHS = UpdateAndGetConstraintRHS(mTime) - bRHS;
        auto extForce = CalculateCurrentExternalLoad(mTime);

        mStructure->GetLogger() << "TimeStep: " << mTimeStep << '\n';


        // extrapolate the history variables
        //            ExtrapolateStaticData(timeStep);


        for (const auto& activeDofSet : mStepActiveDofs)
        {
            mStructure->DofTypeSetIsActive(activeDofSet);


            // if a single active dof is in mDofsWithConstantHessian
            bool usePreFactorizedHessian =
                    activeDofSet.size() == 1 &&
                    mDofsWithConstantHessian.find(*activeDofSet.begin()) != mDofsWithConstantHessian.end();

            calculateStaticData.SetCalculateStaticData(eCalculateStaticData::EULER_FORWARD);
            //                calculateStaticData.SetIndexOfPreviousStaticData(0);

            if (usePreFactorizedHessian)
            {
                auto activeDof = *activeDofSet.begin();

                mStructure->Evaluate(input, evalInternalGradient); // internal gradient only

                Eigen::VectorXd solution;
                preFactorizedHessians[activeDof].Solution(-intForce.J[activeDof], solution);

                dof_dt0.J[activeDof] += solution;
            }
            else
            {
                mStructure->Evaluate(input, evalInternalGradientAndHessian0);

                delta_dof_dt0.J.SetZero();
                delta_dof_dt0.K = deltaBRHS;

                residual = (extForce - intForce) + prevExtForce + extForce;
                residual -= hessian0 * delta_dof_dt0;
                const auto& CMat = mStructure->GetAssembler().GetConstraintMatrix();
                residual.J = Assembler::ApplyCMatrix(residual, CMat);
                hessian0.ApplyCMatrix(CMat);

                delta_dof_dt0.J = mSolver->Solve(hessian0.JJ, residual.J);
                delta_dof_dt0.K = deltaBRHS - CMat * delta_dof_dt0.J;

                dof_dt0 += delta_dof_dt0;
            }
            mStructure->NodeMergeDofValues(dof_dt0);

        } // end for mStepActiveDofs


        if (mCallback && mCallback->Exit(*mStructure))
            return;

        bool acceptSolution = true;

        if (mAutomaticTimeStepping)
            acceptSolution = CheckExtrapolationAndAdjustTimeStep();

        if (acceptSolution)
        {
            // save the new implicit history variables
            calculateStaticData.SetCalculateStaticData(eCalculateStaticData::EULER_BACKWARD);
            calculateStaticData.SetIndexOfPreviousStaticData(1);
            //                mStructure->Evaluate(input, evalUpdateStaticData);

            mStructure->DofTypeActivateAll();
            //                calculateStaticData.SetIndexOfPreviousStaticData(0);
            //                calculateStaticData.SetCalculateStaticData(eCalculateStaticData::USE_PREVIOUS);
            mStructure->Evaluate(input, evalInternalGradient);


            if (mTime + mTimeStep > rTimeDelta)
                mTimeStep = rTimeDelta - mTime;

            mStructure->ElementTotalShiftStaticDataToPast(); // shift static data by one to the past
            timeStep.SetCurrentTimeStep(mTimeStep); // shift time steps by one to the past


            residual = extForce - intForce;
            mPostProcessor->PostProcess(residual);

            lastConverged_dof_dt0 = dof_dt0;

            ++mIterationCount;
        }
        else
        {
            // do not save static data
            // do not shift the time but set the new time step
            mTime -= timeStep[0];
            timeStep[0] = mTimeStep;


            // restore the old solution
            mStructure->DofTypeActivateAll();
            dof_dt0 = lastConverged_dof_dt0;
            mStructure->NodeMergeDofValues(dof_dt0);

            ++mIterationCount;
        }

    } // end while

    for (auto& itPair : preFactorizedHessians)
    {
        itPair.second.CleanUp();
    }
}

double NuTo::ImplicitExplicitBase::CalculateCriticalTimeStep() const
{
    return M_PI;
}

void NuTo::ImplicitExplicitBase::AddDofWithConstantHessian0(Node::eDof rDofWithConstantHessian)
{
    mDofsWithConstantHessian.insert(rDofWithConstantHessian);
}

void NuTo::ImplicitExplicitBase::FactorizeConstantHessians(
        std::map<Node::eDof, SparseDirectSolverMUMPS>& rPreFactorizedHessians)
{
    for (auto dof : mDofsWithConstantHessian)
    {
        if (mStructure->GetNumDependentDofs(dof) > 0)
            throw Exception(__PRETTY_FUNCTION__, "Constant Hessian with constrained dofs is currently not supported.");


        mStructure->DofTypeSetIsActive(dof, true);
        auto hessian0 = mStructure->BuildGlobalHessian0();
        auto hessian0_CSR = hessian0.JJ.ExportToCSRGeneral();
        hessian0_CSR.SetOneBasedIndexing();

        rPreFactorizedHessians[dof].Factorization(hessian0_CSR);
        rPreFactorizedHessians[dof].SetShowTime(false);
        mStructure->DofTypeSetIsActive(dof, false);
    }
}
