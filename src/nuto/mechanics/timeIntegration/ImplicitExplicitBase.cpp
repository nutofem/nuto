/*
 * ImplicitExplicitBase.cpp
 *
 *  Created on: 17 Mar 2016
 *      Author: ttitsche
 */

#include "nuto/mechanics/timeIntegration/ImplicitExplicitBase.h"
#include "nuto/base/Timer.h"
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/structures/StructureOutputDummy.h"

#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveTimeStep.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveCalculateStaticData.h"
#include "nuto/mechanics/constitutive/staticData/ConstitutiveStaticDataGradientDamage.h"


NuTo::ImplicitExplicitBase::ImplicitExplicitBase(StructureBase* rStructure) : TimeIntegrationBase(rStructure)
{}

NuTo::Error::eError NuTo::ImplicitExplicitBase::Solve(double rTimeDelta)
{
    NuTo::Timer timerFull(__PRETTY_FUNCTION__, mStructure->GetShowTime(), mStructure->GetLogger());

    try
    {
        NuTo::Timer timerDebug("Init", true, mStructure->GetLogger());

        mStructure->NodeBuildGlobalDofs(__PRETTY_FUNCTION__);

        CalculateStaticAndTimeDependentExternalLoad();

        const DofStatus& dofStatus = mStructure->GetDofStatus();

        if (mStepActiveDofs.empty())
            throw MechanicsException(__PRETTY_FUNCTION__, "Define a set of active dofs for each calculation step. ");

        // deactivate all dof types
        mStructure->DofTypeDeactivateAll();

        std::map<Node::eDof, SparseDirectSolverMUMPS> preFactorizedHessians;
        FactorizeConstantHessians(preFactorizedHessians);


        /*---------------------------------*\
        |        Allocate Variables         |
        \*---------------------------------*/

        StructureOutputBlockMatrix  hessian0(dofStatus, true);

        StructureOutputBlockVector  delta_dof_dt0(dofStatus, true);

        StructureOutputBlockVector  dof_dt0(dofStatus, true); // e.g. disp
        StructureOutputBlockVector  lastConverged_dof_dt0(dofStatus, true); // e.g. disp

        StructureOutputBlockVector  extForce(dofStatus, true);
        StructureOutputBlockVector  intForce(dofStatus, true);
        StructureOutputBlockVector  residual(dofStatus, true);

        ConstitutiveTimeStep<2> timeStep;
        ConstitutiveCalculateStaticData calculateStaticData(CalculateStaticData::USE_PREVIOUS ,1);

        // for constraints
        // ---------------

        BlockFullVector<double> residual_mod(dofStatus);

        /*---------------------------------*\
        |    Declare and fill Output Maps   |
        \*---------------------------------*/

        // Declare output maps
        std::map<NuTo::StructureEnum::eOutput, NuTo::StructureOutputBase*> evalUpdateStaticData;
        StructureOutputDummy dummy;
        evalUpdateStaticData                [StructureEnum::UPDATE_STATIC_DATA] = &dummy;

        std::map<NuTo::StructureEnum::eOutput, NuTo::StructureOutputBase*> evalInternalGradient;
        evalInternalGradient                [StructureEnum::INTERNAL_GRADIENT] = &intForce;

        std::map<NuTo::StructureEnum::eOutput, NuTo::StructureOutputBase*> evalInternalGradientAndHessian0;
        evalInternalGradientAndHessian0     [StructureEnum::INTERNAL_GRADIENT] = &intForce;
        evalInternalGradientAndHessian0     [StructureEnum::HESSIAN0] = &hessian0;

        /*---------------------------------*\
        |    Declare and fill Input map     |
        \*---------------------------------*/



        ConstitutiveInputMap input;
        input[Constitutive::Input::TIME_STEP] = &timeStep;
        input[Constitutive::Input::CALCULATE_STATIC_DATA] = &calculateStaticData;

        timeStep.SetCurrentTimeStep(mTimeStep);
        timeStep.SetCurrentTimeStep(mTimeStep);
        timeStep.SetCurrentTimeStep(mTimeStep);
        timeStep.SetCurrentTimeStep(mTimeStep);
        timeStep.SetCurrentTimeStep(mTimeStep);

        int gElementsTotal = dynamic_cast<Structure*>(mStructure)->GroupGetElementsTotal();
        auto allElementIds = mStructure->GroupGetMemberIds(gElementsTotal);

        while (mTime < rTimeDelta)
        {
            timerDebug.Reset("current time: " + std::to_string(mTime));

            mStructure->DofTypeSetIsActive(dofStatus.GetDofTypes()); // activate all


            // calculate Delta_BRhs and Delta_ExtForce
            auto bRHS = UpdateAndGetConstraintRHS(mTime);
            auto prevExtForce = CalculateCurrentExternalLoad(mTime);

            mTime += timeStep[0];


            auto deltaBRHS = UpdateAndGetConstraintRHS(mTime) - bRHS;
            auto extForce = CalculateCurrentExternalLoad(mTime);

            timeStep.SetCurrentTimeStep(mTimeStep);

            std::cout << "TimeStep: " << mTimeStep << std::endl;

            // extrapolate the history variables
            calculateStaticData.SetCalculateStaticData(CalculateStaticData::EULER_FORWARD);
            mStructure->Evaluate(input, evalUpdateStaticData);

            for (const auto& activeDofSet : mStepActiveDofs)
            {
                mStructure->DofTypeSetIsActive(activeDofSet);

                // hardcode version for gradient damage model:
                if (activeDofSet.find(Node::NONLOCALEQSTRAIN) != activeDofSet.end())
                {
                    // nonlocal part, use precomputed factorization
                    calculateStaticData.SetCalculateStaticData(CalculateStaticData::USE_PREVIOUS);
                    mStructure->Evaluate(input, evalInternalGradient);

//                    mStructure->GetLogger() << "Initial trial residual:               " << intForce.J.CalculateInfNorm() << "\n";

                    NuTo::FullVector<double,Eigen::Dynamic> solution;
                    preFactorizedHessians[Node::NONLOCALEQSTRAIN].Solution(intForce.J[Node::NONLOCALEQSTRAIN], solution);
                    delta_dof_dt0.SetZero();
                    delta_dof_dt0.J.Import(-solution);


                    dof_dt0 += delta_dof_dt0;

                    mStructure->NodeMergeDofValues(dof_dt0);
                }
                else
                {
                    // displacement part, calculate everything with old data
                    assert(activeDofSet.find(Node::DISPLACEMENTS) != activeDofSet.end());

                    calculateStaticData.SetCalculateStaticData(CalculateStaticData::USE_PREVIOUS);
                    mStructure->Evaluate(input, evalInternalGradientAndHessian0);

                    delta_dof_dt0.J.SetZero();
                    delta_dof_dt0.K = deltaBRHS;

                    residual = intForce;// - prevExtForce - extForce;
                    residual += hessian0 * delta_dof_dt0;
                    residual.ApplyCMatrix(mStructure->GetConstraintMatrix());

//                    mStructure->GetLogger() << "Initial trial residual:               " << residual.J.CalculateInfNorm() << "\n";

                    hessian0.ApplyCMatrix(mStructure->GetConstraintMatrix());


                    delta_dof_dt0.J =  mStructure->SolveBlockSystem(hessian0.JJ, residual.J);
                    delta_dof_dt0.K = deltaBRHS;// - mStructure->GetConstraintMatrix()*delta_dof_dt0.J;



                    dof_dt0 += delta_dof_dt0;

                    mStructure->NodeMergeDofValues(dof_dt0);
//
//                    calculateStaticData.SetCalculateStaticData(CalculateStaticData::EULER_FORWARD);
//                    mStructure->Evaluate(input, evalInternalGradient);
//
//                    std::cout << intForce.J.CalculateInfNorm() << std::endl;
//                    if (intForce.J.CalculateInfNorm()[Node::DISPLACEMENTS] > 1.e-8)
//                    {
//                        std::cout << intForce.J.CalculateInfNorm() << std::endl;
//                        throw MechanicsException("FIX THE SOLVE, YOU IDIOT. I MEAN ITS A LINEAR SYSTEM. HOW HARD CAN IT BE TO SOLVE THAT SHIT PROPERLY.");
//                    }
                }
            } // end for mStepActiveDofs


            // save the new implicit history variables
            calculateStaticData.SetCalculateStaticData(CalculateStaticData::EULER_BACKWARD);
            calculateStaticData.SetIndexOfPreviousStaticData(1);
            mStructure->Evaluate(input, evalUpdateStaticData);


            mStructure->ElementTotalUpdateStaticData();
            mStructure->ElementTotalSaveStaticData();
            PostProcess(residual);

            if (mCallback && mCallback->Exit(*mStructure))
                return Error::SUCCESSFUL;

        } // end while




    }
    catch (MechanicsException& e)
    {
        e.AddMessage(__PRETTY_FUNCTION__, " ERROR performing IMPL-EX.");
        throw e;
    }
    return NuTo::Error::SUCCESSFUL;

}

double NuTo::ImplicitExplicitBase::CalculateCriticalTimeStep() const
{
    return M_PI;
}


void NuTo::ImplicitExplicitBase::AddCalculationStep(const std::set<NuTo::Node::eDof>& rActiveDofs)
{
    mStepActiveDofs.push_back(rActiveDofs);
}

void NuTo::ImplicitExplicitBase::AddDofWithConstantHessian(Node::eDof rDofWithConstantHessian)
{
    mDofsWithConstantHessian.insert(rDofWithConstantHessian);
}

void NuTo::ImplicitExplicitBase::FactorizeConstantHessians(std::map<Node::eDof, SparseDirectSolverMUMPS>& rPreFactorizedHessians)
{
    for (auto dof : mDofsWithConstantHessian)
    {
        if (mStructure->GetNumDependentDofs(dof) > 0)
            throw MechanicsException(__PRETTY_FUNCTION__, "Constant Hessian with constrained dofs is currently not supported.");


        mStructure->DofTypeSetIsActive(dof, true);
        auto hessian0 = mStructure->BuildGlobalHessian0();
        auto hessian0_CSR = hessian0.JJ.ExportToCSRGeneral();
        hessian0_CSR.SetOneBasedIndexing();

        rPreFactorizedHessians[dof].Factorization(hessian0_CSR);
        mStructure->DofTypeSetIsActive(dof, false);
    }
}
