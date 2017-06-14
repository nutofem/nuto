#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp> #include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/ptr_container/serialize_ptr_map.hpp>
#endif // ENABLE_SERIALIZATION

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

//! @brief constructor
//! @param mDimension number of nodes
NuTo::NewmarkDirect::NewmarkDirect(StructureBase* rStructure)
    : NewmarkBase(rStructure)

{
}


//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::NewmarkDirect::Info() const
{
    NewmarkBase::Info();
}

void NuTo::NewmarkDirect::Solve(double rTimeDelta)
{
    NuTo::Timer timerFull(__FUNCTION__, mStructure->GetShowTime(), mStructure->GetLogger());

    // renumber dofs and build constraint matrix
    mStructure->NodeBuildGlobalDofs(__PRETTY_FUNCTION__);

    CalculateStaticAndTimeDependentExternalLoad();


    mToleranceResidual.DefineDefaultValueToIninitializedDofTypes(mToleranceForce);

    if (mMaxTimeStep == 0)
        throw MechanicsException(__PRETTY_FUNCTION__, "max time step is set to zero.");

    double curTime = mTime;
    double timeStep = mTimeStep;

    const DofStatus& dofStatus = mStructure->GetDofStatus();

    /*---------------------------------*\
    | Check number of calculation steps |
    \*---------------------------------*/

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


    /*---------------------------------*\
    |        Allocate Variables         |
    \*---------------------------------*/

    StructureOutputBlockMatrix hessian0(dofStatus, true);
    StructureOutputBlockMatrix hessian1(dofStatus);
    StructureOutputBlockMatrix hessian2(dofStatus);

    // Vectors
    // -------

    StructureOutputBlockVector delta_dof_dt0(dofStatus, true);

    StructureOutputBlockVector dof_dt0(dofStatus, true); // e.g. disp
    StructureOutputBlockVector dof_dt1(dofStatus, true); // e.g. velocity
    StructureOutputBlockVector dof_dt2(dofStatus, true); // e.g. accelerations

    StructureOutputBlockVector trial_dof_dt0(dofStatus, true); // e.g. disp
    StructureOutputBlockVector trial_dof_dt1(dofStatus, true); // e.g. velocity
    StructureOutputBlockVector trial_dof_dt2(dofStatus, true); // e.g. accelerations

    StructureOutputBlockVector lastConverged_dof_dt0(dofStatus, true); // e.g. disp
    StructureOutputBlockVector lastConverged_dof_dt1(dofStatus, true); // e.g. velocity
    StructureOutputBlockVector lastConverged_dof_dt2(dofStatus, true); // e.g. accelerations


    StructureOutputBlockVector extForce(dofStatus, true);
    StructureOutputBlockVector intForce(dofStatus, true);
    StructureOutputBlockVector residual(dofStatus, true);

    StructureOutputBlockVector prevIntForce(dofStatus, true);
    StructureOutputBlockVector prevExtForce(dofStatus, true);
    StructureOutputBlockVector prevResidual(dofStatus, true);


    // for constraints
    // ---------------

    BlockFullVector<double> residual_mod(dofStatus);
    BlockFullVector<double> bRHS(dofStatus);
    BlockFullVector<double> deltaBRHS(dofStatus);

    const auto& cmat = mStructure->GetAssembler().GetConstraintMatrix();

    /*---------------------------------*\
    |    Declare and fill Output Maps   |
    \*---------------------------------*/

    // Declare output maps
    std::map<NuTo::eStructureOutput, NuTo::StructureOutputBase*> evaluate_InternalGradient;
    std::map<NuTo::eStructureOutput, NuTo::StructureOutputBase*> evaluate_InternalGradient_Hessian0Hessian1;
    std::map<NuTo::eStructureOutput, NuTo::StructureOutputBase*> evaluate_Hessian0_Hessian1;

    evaluate_InternalGradient[eStructureOutput::INTERNAL_GRADIENT] = &intForce;
    evaluate_InternalGradient_Hessian0Hessian1[eStructureOutput::INTERNAL_GRADIENT] = &intForce;

    evaluate_InternalGradient_Hessian0Hessian1[eStructureOutput::HESSIAN0] = &hessian0;
    evaluate_Hessian0_Hessian1[eStructureOutput::HESSIAN0] = &hessian0;

    if (mStructure->GetNumTimeDerivatives() >= 1 && mMuDampingMass == 0.)
    {
        hessian1.Resize(dofStatus.GetNumActiveDofsMap(), dofStatus.GetNumDependentDofsMap());
        evaluate_InternalGradient_Hessian0Hessian1[eStructureOutput::HESSIAN1] = &hessian1;
        evaluate_Hessian0_Hessian1[eStructureOutput::HESSIAN1] = &hessian1;
    }

    if (mStructure->GetNumTimeDerivatives() >= 2)
    {
        hessian2.Resize(dofStatus.GetNumActiveDofsMap(), dofStatus.GetNumDependentDofsMap());
        if (mUseLumpedMass)
        {
            hessian2 = mStructure->BuildGlobalHessian2Lumped();
        }
        else
        {
            hessian2 = mStructure->BuildGlobalHessian2();
        }
        CalculateMuDampingMatrix(hessian1, hessian2);
    }


    /*---------------------------------*\
    |    Declare and fill Input Map     |
    \*---------------------------------*/

    ConstitutiveInputMap inputMap;
    inputMap[Constitutive::eInput::CALCULATE_STATIC_DATA] =
            std::make_unique<ConstitutiveCalculateStaticData>(eCalculateStaticData::EULER_BACKWARD);
    double& inputTime =
            (*inputMap.emplace(Constitutive::eInput::TIME, std::make_unique<ConstitutiveScalar>()).first->second)[0];
    inputTime = mTime;

    ExtractDofValues(lastConverged_dof_dt0, lastConverged_dof_dt1, lastConverged_dof_dt2);

    UpdateAndGetAndMergeConstraintRHS(curTime, lastConverged_dof_dt0);

    // ******************************************************
    mStructure->Evaluate(inputMap, evaluate_InternalGradient_Hessian0Hessian1);
    // ******************************************************

    StructureOutputBlockVector initialExtForce =
            CalculateCurrentExternalLoad(curTime); // put this in element evaluate soon!

    // set first time derivative for temperature problem automatically
    // for (const auto& activeDofs : mStepActiveDofs)
    //{
    //    auto temp_iterator = activeDofs.find(NuTo::Node::eDof::TEMPERATURE);
    //    bool is_temperature = temp_iterator != activeDofs.end();
    //    if (mStructure->GetNumTimeDerivatives() == 1 && is_temperature)
    //    {
    //        auto rhs = hessian0*lastConverged_dof_dt0 - initialExtForce;
    //        lastConverged_dof_dt1.J = mStructure->SolveBlockSystem(hessian1.JJ, rhs.J);
    //        mStructure->NodeMergeDofValues(1, lastConverged_dof_dt1);
    //        mStructure->Evaluate(inputMap, evaluate_InternalGradient_Hessian0Hessian1);
    //    }
    //}

    residual = CalculateResidual(intForce, initialExtForce, hessian2, lastConverged_dof_dt1, lastConverged_dof_dt2);
    residual.ApplyCMatrix(residual_mod, cmat);
    if (mToleranceResidual < residual_mod.CalculateInfNorm())
    {
        mStructure->GetLogger() << residual_mod.CalculateInfNorm();
        throw MechanicsException(__PRETTY_FUNCTION__, "Initial configuration is not in (dynamic) equilibrium.");
    }

    CalculateResidualKForPostprocessing(residual, hessian2, lastConverged_dof_dt1, lastConverged_dof_dt2);
    PostProcess(residual);


    // the minimal time step defined, which is equivalent to six cut-backs
    if (mAutomaticTimeStepping)
    {
        SetMinTimeStep(mMinTimeStep > 0. ? mMinTimeStep : mTimeStep * std::pow(0.5, 6.));
    }

    while (curTime < rTimeDelta)
    {


        //! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //! LEAVE IT OR FIX IT IF YOU CAN:
        //! If you dont make a copy of the dof set with dofStatus.GetDofTypes() and use it directly in the for loop
        //! everything seems
        //! to work fine, but valgrind tells you otherwise. Problem is, that the DofType is changed during the function
        //! call inside
        //! the loop which leads to reads in already freed blocks of memory
        //! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        std::set<Node::eDof> currentDofSet = dofStatus.GetDofTypes();
        for (const auto& dof : currentDofSet)
        {
            mStructure->DofTypeSetIsActive(dof, true);
        }


        if (timeStep < mMinTimeStep)
            throw MechanicsException(__PRETTY_FUNCTION__,
                                     "time step is smaller than minimum - no convergence is obtained.");

        // calculate Delta_BRhs and Delta_ExtForce
        bRHS = UpdateAndGetAndMergeConstraintRHS(curTime, lastConverged_dof_dt0);
        prevExtForce = CalculateCurrentExternalLoad(curTime);

        curTime += timeStep;
        inputTime = mTime + timeStep;
        SetTimeAndTimeStep(
                curTime, timeStep,
                rTimeDelta); // check whether harmonic excitation, check whether curTime is too close to the time data

        deltaBRHS = UpdateAndGetConstraintRHS(curTime) - bRHS;

        unsigned int staggeredStepNumber = 0; // at the moment needed to do the postprocessing after the last step and
        // not after every step of a staggered solution.
        for (const auto& activeDofs : mStepActiveDofs)
        {

            ++staggeredStepNumber;
            mStructure->DofTypeSetIsActive(activeDofs);


            PrintInfoStagger();

            // ******************************************************
            mStructure->Evaluate(inputMap, evaluate_Hessian0_Hessian1);
            // ******************************************************

            extForce = CalculateCurrentExternalLoad(curTime);

            /*------------------------------------------------*\
            |         Calculate Residual for trail state       |
            \*------------------------------------------------*/
            residual = extForce - prevExtForce;
            CalculateResidualTrial(residual, deltaBRHS, hessian0, hessian1, hessian2, lastConverged_dof_dt1,
                                   lastConverged_dof_dt2, timeStep);
            residual.ApplyCMatrix(residual_mod, cmat);

            mStructure->GetLogger() << "\n"
                                    << "Initial trial residual:               " << residual_mod.CalculateInfNorm()
                                    << "\n";

            // ******************************************************
            delta_dof_dt0.J = BuildHessianModAndSolveSystem(hessian0, hessian1, hessian2, residual_mod, timeStep);
            delta_dof_dt0.K = deltaBRHS - cmat * delta_dof_dt0.J;
            ++mIterationCount;
            // ******************************************************

            // calculate trial state
            dof_dt0 = lastConverged_dof_dt0 + delta_dof_dt0;
            if (mStructure->GetNumTimeDerivatives() >= 1)
                dof_dt1 = CalculateDof1(delta_dof_dt0, lastConverged_dof_dt1, lastConverged_dof_dt2, timeStep);
            if (mStructure->GetNumTimeDerivatives() >= 2)
                dof_dt2 = CalculateDof2(delta_dof_dt0, lastConverged_dof_dt1, lastConverged_dof_dt2, timeStep);


            MergeDofValues(dof_dt0, dof_dt1, dof_dt2, false);


            // ******************************************************
            mStructure->Evaluate(inputMap, evaluate_InternalGradient);
            // ******************************************************

            residual = CalculateResidual(intForce, extForce, hessian2, dof_dt1, dof_dt2);
            residual.ApplyCMatrix(residual_mod, cmat);


            BlockScalar normResidual = residual_mod.CalculateInfNorm();

            int iteration = 0;
            while (!(normResidual < mToleranceResidual) && iteration < mMaxNumIterations)
            {

                // ******************************************************
                mStructure->Evaluate(inputMap, evaluate_Hessian0_Hessian1);
                // ******************************************************

                if (mCheckCoefficientMatrix)
                    mStructure->ElementCheckHessian0(1.e-6, 1.e-8);


                // ******************************************************
                delta_dof_dt0.J = BuildHessianModAndSolveSystem(hessian0, hessian1, hessian2, residual_mod, timeStep);
                delta_dof_dt0.K = cmat * delta_dof_dt0.J * (-1.);
                ++mIterationCount;
                // ******************************************************


                // perform a line search
                double alpha = 1;
                BlockScalar trialNormResidual(mStructure->GetDofStatus());
                do
                {
                    // calculate line search trial state
                    trial_dof_dt0 = dof_dt0 + delta_dof_dt0 * alpha;
                    if (mStructure->GetNumTimeDerivatives() >= 1)
                        trial_dof_dt1 = dof_dt1 + delta_dof_dt0 * (alpha * mGamma / (timeStep * mBeta));
                    if (mStructure->GetNumTimeDerivatives() >= 2)
                        trial_dof_dt2 = dof_dt2 + delta_dof_dt0 * (alpha / (timeStep * timeStep * mBeta));

                    MergeDofValues(trial_dof_dt0, trial_dof_dt1, trial_dof_dt2, false);

                    // ******************************************************
                    mStructure->Evaluate(inputMap, evaluate_InternalGradient);
                    // ******************************************************

                    residual = CalculateResidual(intForce, extForce, hessian2, trial_dof_dt1, trial_dof_dt2);
                    residual.ApplyCMatrix(residual_mod, cmat);

                    trialNormResidual = residual_mod.CalculateInfNorm();

                    mStructure->GetLogger() << "[Linesearch a=" << std::to_string(alpha).substr(0, 6)
                                            << "] Trial residual: " << trialNormResidual << "\n";

                    alpha *= 0.5;

                } while (mPerformLineSearch && alpha > mMinLineSearchStep &&
                         trialNormResidual > (1. - alpha) * normResidual);

                if (alpha > mMinLineSearchStep || !mPerformLineSearch)
                {
                    // improvement is achieved, go to next Newton step
                    dof_dt0 = trial_dof_dt0;
                    if (mStructure->GetNumTimeDerivatives() >= 1)
                        dof_dt1 = trial_dof_dt1;
                    if (mStructure->GetNumTimeDerivatives() >= 2)
                        dof_dt2 = trial_dof_dt2;
                    normResidual = trialNormResidual;

                    PrintInfoIteration(normResidual, iteration);
                    iteration++;
                }
                else
                {
                    // and leave
                    iteration = mMaxNumIterations;
                }

            } // end of while(normResidual<mToleranceForce && iteration<mMaxNumIterations)

            if (normResidual < mToleranceResidual)
            {
                // converged solution
                if (mVerboseLevel > 2)
                    mStructure->GetLogger() << "\n *** UpdateStaticData *** from NewMarkDirect \n";

                // Update static data
                mStructure->ElementTotalUpdateStaticData();

                // store converged step
                lastConverged_dof_dt0 = dof_dt0;
                if (mStructure->GetNumTimeDerivatives() >= 1)
                    lastConverged_dof_dt1 = dof_dt1;
                if (mStructure->GetNumTimeDerivatives() >= 1)
                    lastConverged_dof_dt2 = dof_dt2;

                prevResidual = residual;


                MergeDofValues(dof_dt0, dof_dt1, dof_dt2, true);


                mTime += timeStep;
                inputTime = mTime;

                mStructure->GetLogger() << "Convergence after " << iteration << " iterations at time " << mTime
                                        << " (timestep " << timeStep << ").\n";
                mStructure->GetLogger() << "Residual: \t" << normResidual << "\n";
                // perform Postprocessing
                if (staggeredStepNumber >= mStepActiveDofs.size())
                {
                    CalculateResidualKForPostprocessing(prevResidual, hessian2, dof_dt1, dof_dt2);
                    PostProcess(prevResidual);
                }


                // eventually increase next time step
                if (mAutomaticTimeStepping && iteration < 0.25 * mMaxNumIterations)
                {
                    timeStep *= 1.5;
                    if (timeStep > mMaxTimeStep)
                        timeStep = mMaxTimeStep;
                }

                if (mCallback && mCallback->Exit(*mStructure))
                    return;
            }
            else
            {
                mStructure->GetLogger() << "No convergence with timestep " << timeStep << " at time " << mTime << "\n";
                // no convergence
                if (mAutomaticTimeStepping)
                {
                    // no convergence, reduce the time step and start from scratch
                    curTime -= timeStep;
                    timeStep *= 0.5;
                    if (timeStep < mMinTimeStep)
                    {
                        mStructure->GetLogger() << "The minimal time step achieved, the actual time step is "
                                                << timeStep << "\n";
                        throw MechanicsException(__PRETTY_FUNCTION__,
                                                 "No convergence, the current time step is too short.");
                    }
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
    } // while time steps
}


void NuTo::NewmarkDirect::CalculateMuDampingMatrix(StructureOutputBlockMatrix& rHessian_dt1,
                                                   const StructureOutputBlockMatrix& rHessian_dt2) const
{
    if (mMuDampingMass != 0)
    {
        if (mStructure->GetNumTimeDerivatives() < 2)
            throw MechanicsException(std::string("[") + __PRETTY_FUNCTION__ +
                                     "] MuDampingMass requires a mass matrix (2nd time derivatives).");
        //        if (!rHessian2.IsConstant())
        //            throw MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] MuDampingMass requires a
        //            constant mass matrix.");

        rHessian_dt1.Resize(mStructure->GetDofStatus().GetNumActiveDofsMap(),
                            mStructure->GetDofStatus().GetNumDependentDofsMap());
        rHessian_dt1.SetZero();
        rHessian_dt1.AddScal(rHessian_dt2, mMuDampingMass);
    }
}

NuTo::StructureOutputBlockVector NuTo::NewmarkDirect::CalculateDof1(const StructureOutputBlockVector& rDeltaDof_dt0,
                                                                    const StructureOutputBlockVector& rDof_dt1,
                                                                    const StructureOutputBlockVector& rDof_dt2,
                                                                    double rTimeStep) const
{
    return rDeltaDof_dt0 * (mGamma / (rTimeStep * mBeta)) + rDof_dt1 * (1. - mGamma / mBeta) +
           rDof_dt2 * (rTimeStep * (1. - mGamma / (2. * mBeta)));
}


NuTo::StructureOutputBlockVector NuTo::NewmarkDirect::CalculateDof2(const StructureOutputBlockVector& rDeltaDof_dt0,
                                                                    const StructureOutputBlockVector& rDof_dt1,
                                                                    const StructureOutputBlockVector& rDof_dt2,
                                                                    double rTimeStep) const
{
    return rDeltaDof_dt0 * (1. / (rTimeStep * rTimeStep * mBeta)) - rDof_dt1 * (1. / (rTimeStep * mBeta)) -
           rDof_dt2 * ((0.5 - mBeta) / mBeta);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

NuTo::StructureOutputBlockVector NuTo::NewmarkDirect::CalculateResidual(
        const StructureOutputBlockVector& rIntForce, const StructureOutputBlockVector& rExtForce,
        const StructureOutputBlockMatrix& rHessian_dt2, const StructureOutputBlockVector& rDof_dt1,
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NuTo::NewmarkDirect::CalculateResidualKForPostprocessing(StructureOutputBlockVector& rResidual,
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NuTo::NewmarkDirect::CalculateResidualTrial(StructureOutputBlockVector& rResidual,
                                                 const BlockFullVector<double>& rDeltaBRHS,
                                                 const StructureOutputBlockMatrix& rHessian_dt0,
                                                 const StructureOutputBlockMatrix& rHessian_dt1,
                                                 const StructureOutputBlockMatrix& rHessian_dt2,
                                                 const StructureOutputBlockVector& rDof_dt1,
                                                 const StructureOutputBlockVector& rDof_dt2, double rTimeStep) const
{


    StructureOutputBlockVector deltaDof_dt0(mStructure->GetDofStatus(), true);
    deltaDof_dt0.J.SetZero();
    deltaDof_dt0.K = rDeltaBRHS;


    rResidual -= rHessian_dt0 * deltaDof_dt0;

    if (mStructure->GetNumTimeDerivatives() >= 1)
    {
        StructureOutputBlockVector delta_dof1 = CalculateDof1(deltaDof_dt0, rDof_dt1, rDof_dt2, rTimeStep) - rDof_dt1;
        rResidual -= rHessian_dt1 * delta_dof1;
    }

    if (mStructure->GetNumTimeDerivatives() >= 2)
    {
        StructureOutputBlockVector delta_dof2 = CalculateDof2(deltaDof_dt0, rDof_dt1, rDof_dt2, rTimeStep) - rDof_dt2;
        rResidual -= rHessian_dt2 * delta_dof2;
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NuTo::NewmarkDirect::PrintInfoStagger() const
{
    if (mStepActiveDofs.size() == 1) // equals unstaggered solution, no info needed
        return;

    NuTo::Logger& logger = mStructure->GetLogger();
    logger << "\n"
           << "Activated Dofs: | ";
    for (auto dof : mStructure->GetDofStatus().GetActiveDofTypes())
    {
        logger << Node::DofToString(dof) << " | ";
    }
    logger << "\n";
}


//! @brief Prints Info about the current iteration
void NuTo::NewmarkDirect::PrintInfoIteration(const BlockScalar& rNormResidual, int rIteration) const
{
    if (mStructure->GetDofStatus().GetDofTypes().find(Node::eDof::NONLOCALEQSTRAIN) !=
        mStructure->GetDofStatus().GetDofTypes().end())
        return; // Hi Volker. Nothing wrong with the 1 line iteration output in the Solve() method. IMO.
    // And please consider line search...

    NuTo::Logger& logger = mStructure->GetLogger();

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


NuTo::BlockFullVector<double>
NuTo::NewmarkDirect::BuildHessianModAndSolveSystem(StructureOutputBlockMatrix& rHessian_dt0,
                                                   const StructureOutputBlockMatrix& rHessian_dt1,
                                                   const StructureOutputBlockMatrix& rHessian_dt2,
                                                   const BlockFullVector<double>& rResidualMod, double rTimeStep) const
{
    Timer timer(__FUNCTION__, GetShowTime(), mStructure->GetLogger());

    //    if (rHessian_dt0.IsConstant())
    //    {
    //        // save the values in rHessian0 for later evaluation and allocate a new hessian object
    //        StructureOutputBlockMatrix hessian(mStructure->GetDofStatus());
    //        hessian = rHessian_dt0;
    //        if (mStructure->GetNumTimeDerivatives() >= 1)
    //            hessian.AddScal(rHessian_dt1, mGamma / (mBeta * rTimeStep));
    //
    //        if (mStructure->GetNumTimeDerivatives() >= 2)
    //            hessian.AddScal(rHessian_dt2, 1. / (mBeta * rTimeStep * rTimeStep));
    //
    //        hessian.ApplyCMatrix(mStructure->GetConstraintMatrix());
    //        return mStructure->SolveBlockSystem(hessian.JJ, rResidualMod);
    //    }
    //    else
    {
        // since rHessian0 will change in the next iteration, the rHessian0 will be the hessian for the solver
        if (mStructure->GetNumTimeDerivatives() >= 1)
            rHessian_dt0.AddScal(rHessian_dt1, mGamma / (mBeta * rTimeStep));

        if (mStructure->GetNumTimeDerivatives() >= 2)
            rHessian_dt0.AddScal(rHessian_dt2, 1. / (mBeta * rTimeStep * rTimeStep));

        rHessian_dt0.ApplyCMatrix(mStructure->GetAssembler().GetConstraintMatrix());

        auto result = mSolver->Solve(rHessian_dt0.JJ, rResidualMod);
        return result;
    }
}


#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
// serializes the class
template void NuTo::NewmarkDirect::serialize(boost::archive::binary_oarchive& ar, const unsigned int version);
template void NuTo::NewmarkDirect::serialize(boost::archive::xml_oarchive& ar, const unsigned int version);
template void NuTo::NewmarkDirect::serialize(boost::archive::text_oarchive& ar, const unsigned int version);
template void NuTo::NewmarkDirect::serialize(boost::archive::binary_iarchive& ar, const unsigned int version);
template void NuTo::NewmarkDirect::serialize(boost::archive::xml_iarchive& ar, const unsigned int version);
template void NuTo::NewmarkDirect::serialize(boost::archive::text_iarchive& ar, const unsigned int version);
template <class Archive>
void NuTo::NewmarkDirect::serialize(Archive& ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialization of NewmarkDirect"
              << "\n";
#endif
    ar& BOOST_SERIALIZATION_BASE_OBJECT_NVP(NewmarkBase) & BOOST_SERIALIZATION_NVP(mMinLineSearchStep);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialization of NewmarkDirect"
              << "\n";
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::NewmarkDirect)
#endif // SWIG
#endif // ENABLE_SERIALIZATION
