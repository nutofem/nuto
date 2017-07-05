#include "math/SparseMatrixCSRVector2.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/groups/Group.h"
#include "mechanics/structures/StructureBase.h"
#include "mechanics/structures/StructureOutputBlockMatrix.h"
#include "mechanics/timeIntegration/VelocityVerlet.h"
#include "mechanics/timeIntegration/TimeIntegrationEnum.h"
#include "mechanics/structures/Assembler.h"

#include "base/Timer.h"

//! @brief constructor
//! @param mDimension number of nodes
NuTo::VelocityVerlet::VelocityVerlet(StructureBase* rStructure)
    : TimeIntegrationBase(rStructure)
{
}


//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::VelocityVerlet::Info() const
{
    TimeIntegrationBase::Info();
}

//! @brief calculate the critical time step for explicit routines
//! for implicit routines, this will simply return zero (cmp HasCriticalTimeStep())
double NuTo::VelocityVerlet::CalculateCriticalTimeStep() const
{
    double maxGlobalEigenValue = mStructure->ElementTotalCalculateLargestElementEigenvalue();
    return 2. / std::sqrt(maxGlobalEigenValue);
}


//! @brief perform the time integration
//! @param rStructure ... structure
//! @param rTimeDelta ... length of the simulation
void NuTo::VelocityVerlet::Solve(double rTimeDelta)
{
    NuTo::Timer timer(__FUNCTION__, mStructure->GetShowTime(), mStructure->GetLogger());

    try
    {
        mStructure->NodeBuildGlobalDofs(__PRETTY_FUNCTION__);

        if (mStructure->GetDofStatus().HasInteractingConstraints())
            throw MechanicsException(__PRETTY_FUNCTION__,
                                     "not implemented for constrained systems including multiple dofs.");

        if (mTimeStep == 0.)
        {
            if (this->HasCriticalTimeStep())
            {
                mTimeStep = this->CalculateCriticalTimeStep();
            }
            else
            {
                throw MechanicsException(
                        "[NuTo::VelocityVerlet::Solve] time step not set for unconditional stable algorithm.");
            }
        }

        // check constraints:
        mStructure->GetAssembler().ConstraintUpdateRhs(42);
        if (mStructure->GetAssembler().GetConstraintRhs().Export().sum() > 1.e-6)
        {
            throw MechanicsException(__PRETTY_FUNCTION__, "solution with constraints != 0 not yet implemented.");
        }
        mStructure->GetAssembler().ConstraintUpdateRhs(0);

        std::cout << "time step " << mTimeStep << std::endl;
        std::cout << "number of time steps " << rTimeDelta / mTimeStep << std::endl;


        StructureOutputBlockVector outOfBalance(mStructure->GetDofStatus(), true);

        CalculateStaticAndTimeDependentExternalLoad();

        // store last converged displacements, velocities and accelerations
        auto dof_dt0 = mStructure->NodeExtractDofValues(0);
        auto dof_dt1 = mStructure->NodeExtractDofValues(1);
        auto dof_dt2 = mStructure->NodeExtractDofValues(2);

        // calculate lumped mass matrix
        StructureOutputBlockMatrix hessian2 = mStructure->BuildGlobalHessian2Lumped();
        double numericMass = 0;
        numericMass += hessian2.JJ(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS).Sum();
        numericMass += hessian2.JK(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS).Sum();
        numericMass += hessian2.KJ(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS).Sum();
        numericMass += hessian2.KK(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS).Sum();

        numericMass /= mStructure->GetDimension(); // since the mass is added to nodes in every direction
        std::cout << "the total mass is " << numericMass << std::endl;

        // invert the mass matrix
        hessian2.CwiseInvert();


        double curTime = 0.;

        auto extLoad = CalculateCurrentExternalLoad(curTime);
        auto intForce = mStructure->BuildGlobalInternalGradient();

        while (curTime < rTimeDelta)
        {
            // increase time step
            curTime += mTimeStep;
            if (mStructure->GetVerboseLevel() > 5)
                std::cout << "curTime " << curTime << " (" << curTime / rTimeDelta
                          << ") max Disp = " << dof_dt0.J[Node::eDof::DISPLACEMENTS].maxCoeff() << std::endl;

            // apply constraints for new time stepdouble RHSConstraint
            //            double timeDependentConstraintFactor(0);
            // calculate new displacement approximation
            dof_dt0 += dof_dt1 * mTimeStep + dof_dt2 * (mTimeStep * mTimeStep * 0.5);

            dof_dt0.K = mStructure->NodeCalculateDependentDofValues(dof_dt0.J); //?
            mStructure->NodeMergeDofValues(0, dof_dt0);
            if (mMergeActiveDofValuesOrder1)
            {
                dof_dt1.K = mStructure->NodeCalculateDependentDofValues(dof_dt1.J); //?
                mStructure->NodeMergeDofValues(1, dof_dt1);
            }
            if (mMergeActiveDofValuesOrder2)
            {
                dof_dt2.K = mStructure->NodeCalculateDependentDofValues(dof_dt2.J); //?
                mStructure->NodeMergeDofValues(2, dof_dt2);
            }
            mStructure->ElementTotalUpdateTmpStaticData();
            mStructure->ElementTotalUpdateStaticData();

            // calculate external force
            extLoad = CalculateCurrentExternalLoad(curTime);

            // calculate internal force (with update of history variables = true)
            intForce = mStructure->BuildGlobalInternalGradient();

            //**********************************************
            // PostProcessing
            //**********************************************

            // postprocess data for plotting
            this->PostProcess(extLoad - intForce);

            // calculate new accelerations and velocities of independent dofs
            auto dof_dt2_new = dof_dt2;
            dof_dt2_new = hessian2 * (extLoad - intForce);

            dof_dt1 += (dof_dt2 + dof_dt2_new) * (0.5 * mTimeStep);

            // update acceleration
            dof_dt2 = dof_dt2_new;

            // update time, this is different from curTime if several load cycles are sequentially added
            mTime += mTimeStep;
        }
    }
    catch (MechanicsException& e)
    {
        e.AddMessage("[NuTo::VelocityVerlet::Solve] performing Newton-Raphson iteration.");
        throw;
    }
}
