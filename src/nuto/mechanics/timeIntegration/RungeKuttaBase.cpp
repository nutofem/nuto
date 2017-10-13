// $Id: RungeKuttaBase.cpp 575 2011-09-20 18:05:35Z unger3 $

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/ptr_container/serialize_ptr_map.hpp>
#endif // ENABLE_SERIALIZATION

#ifdef _OPENMP
#include <omp.h>
#endif

#include "nuto/math/SparseMatrixCSRVector2.h"
#include "nuto/base/ErrorEnum.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/nodes/NodeEnum.h"
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/structures/StructureOutputBlockMatrix.h"
#include "nuto/mechanics/timeIntegration/RungeKuttaBase.h"
#include "nuto/mechanics/timeIntegration/TimeIntegrationEnum.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/base/Timer.h"

#include "nuto/math/SparseDirectSolverPardiso.h"
#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/math/SparseMatrixCSRSymmetric.h"
#include "nuto/math/SparseMatrixCSRVector2General.h"

//! @brief constructor
//! @param mDimension number of nodes
NuTo::RungeKuttaBase::RungeKuttaBase(StructureBase* rStructure)
    : TimeIntegrationBase(rStructure)
{
    mTimeStep = 0.;
}


//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::RungeKuttaBase::Info() const
{
    TimeIntegrationBase::Info();
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::RungeKuttaBase::serialize(boost::archive::binary_oarchive& ar, const unsigned int version);
template void NuTo::RungeKuttaBase::serialize(boost::archive::xml_oarchive& ar, const unsigned int version);
template void NuTo::RungeKuttaBase::serialize(boost::archive::text_oarchive& ar, const unsigned int version);
template void NuTo::RungeKuttaBase::serialize(boost::archive::binary_iarchive& ar, const unsigned int version);
template void NuTo::RungeKuttaBase::serialize(boost::archive::xml_iarchive& ar, const unsigned int version);
template void NuTo::RungeKuttaBase::serialize(boost::archive::text_iarchive& ar, const unsigned int version);
template <class Archive>
void NuTo::RungeKuttaBase::serialize(Archive& ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialization of RungeKuttaBase"
              << "\n";
#endif
    ar& BOOST_SERIALIZATION_BASE_OBJECT_NVP(TimeIntegrationBase) & BOOST_SERIALIZATION_NVP(mTimeStep);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialization of RungeKuttaBase"
              << "\n";
#endif
}

#endif // ENABLE_SERIALIZATION


//! @brief perform the time integration
//! @param rStructure ... structure
//! @param rTimeDelta ... length of the simulation
NuTo::eError NuTo::RungeKuttaBase::Solve(double rTimeDelta)
{
    NuTo::Timer timer(__PRETTY_FUNCTION__, mStructure->GetShowTime(), mStructure->GetLogger());

    try
    {
        mStructure->NodeBuildGlobalDofs(__PRETTY_FUNCTION__);

        //        if (mStructure->GetDofStatus().HasInteractingConstraints())
        //            throw MechanicsException(__PRETTY_FUNCTION__, "not implemented for constrained systems including
        //            multiple dofs.");

        if (mTimeStep == 0.)
        {
            if (this->HasCriticalTimeStep())
            {
                mTimeStep = this->CalculateCriticalTimeStep();
            }
            else
            {
                throw MechanicsException(
                        "[NuTo::RungeKuttaBase::Solve] time step not set for unconditional stable algorithm.");
            }
        }

        std::cout
                << "modify computation of critical time step to include the dependence on the time integration scheme."
                << std::endl;
        // calculate instead the smallest eigenfrequency, depending on the time integration this gives the critical time
        // step

        std::cout << "time step " << mTimeStep << std::endl;
        std::cout << "number of time steps " << rTimeDelta / mTimeStep << std::endl;

        StructureOutputBlockVector outOfBalance(mStructure->GetDofStatus(), true);

        CalculateStaticAndTimeDependentExternalLoad();

        // store last converged displacements, velocities and accelerations
        auto dof_dt0 = mStructure->NodeExtractDofValues(0);
        auto dof_dt1 = mStructure->NodeExtractDofValues(1);

        StructureOutputBlockMatrix hessian2 = mStructure->BuildGlobalHessian2Lumped();
        double numericMass = 0;
        numericMass += hessian2.JJ(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS).Sum();
        numericMass += hessian2.JK(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS).Sum();
        numericMass += hessian2.KJ(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS).Sum();
        numericMass += hessian2.KK(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS).Sum();

        numericMass /= mStructure->GetDimension(); // since the mass is added to nodes in every direction
        std::cout << "the total mass is " << numericMass << std::endl;

        // invert the mass matrix ==> only if M is diagonal, better solve system (if diagonal, then no additional
        // effort)
        //        hessian2.CwiseInvert();

        hessian2.ApplyCMatrix(mStructure->GetConstraintMatrix());

        double curTime = 0;
        auto extLoad = CalculateCurrentExternalLoad(curTime);
        auto intForce = mStructure->BuildGlobalInternalGradient();

        std::vector<StructureOutputBlockVector> d_dof_dt0_tmp(this->GetNumStages(), mStructure->GetDofStatus());
        std::vector<StructureOutputBlockVector> d_dof_dt1_tmp(this->GetNumStages(), mStructure->GetDofStatus());

// allocate solver
#if defined(HAVE_PARDISO) && defined(_OPENMP)
        NuTo::SparseDirectSolverPardiso mySolver(this->mStructure->GetNumProcessors(),
                                                 this->mStructure->GetVerboseLevel()); // note: not the MKL version
        const BlockSparseMatrix& rMatrix = hessian2.JJ;
        NuTo::SparseMatrixCSRGeneral<double> matrixForSolver = rMatrix.ExportToCSRGeneral();
        matrixForSolver.SetOneBasedIndexing();
        mySolver.Factorization(matrixForSolver);
#else
        NuTo::SparseDirectSolverMUMPS mySolver;
        const BlockSparseMatrix& rMatrix = hessian2.JJ;
        NuTo::SparseMatrixCSRGeneral<double> matrixForSolver = rMatrix.ExportToCSRGeneral();
        matrixForSolver.SetOneBasedIndexing();
        mySolver.Factorization(matrixForSolver);
#endif
        std::vector<double> stageDerivativeFactor(this->GetNumStages() - 1);
        while (curTime < rTimeDelta)
        {
            // calculate for delta_t = 0
            std::cout << "curTime " << curTime << " (" << curTime / rTimeDelta
                      << ") max Disp = " << dof_dt0.J[Node::eDof::DISPLACEMENTS].maxCoeff() << std::endl;
            auto dof_dt0_new = dof_dt0;
            auto dof_dt1_new = dof_dt1;

            double prevTime(mTime);
            double prevCurTime(curTime);
            for (int countStage = 0; countStage < this->GetNumStages(); countStage++)
            {
                // std::cout << "\n stage weight " << GetStageWeights(countStage) << std::endl;
                double deltaTimeStage = this->GetStageTimeFactor(countStage) * mTimeStep;
                this->GetStageDerivativeFactor(stageDerivativeFactor, countStage);
                auto dof_dt0_tmp = dof_dt0;
                auto dof_dt1_tmp = dof_dt1;
                for (int countStage2 = 0; countStage2 < countStage; countStage2++)
                {
                    if (stageDerivativeFactor[countStage2] != 0.)
                    {
                        dof_dt0_tmp.J += d_dof_dt0_tmp[countStage2].J * (stageDerivativeFactor[countStage2]);
                        dof_dt1_tmp.J += d_dof_dt1_tmp[countStage2].J * (stageDerivativeFactor[countStage2]);
                    }
                }

                if (this->HasTimeChanged(countStage) == true)
                {
                    curTime = prevCurTime + deltaTimeStage;
                    mTime = prevTime + deltaTimeStage;

                    // to be implemented mStructure->SetCurrentTime(mTime);
                    // an update of the external load factor and the time dependent constraint is only
                    // necessary for a modified global time
                    if (mTimeDependentConstraint != -1)
                    {
                        throw MechanicsException(
                                "[NuTo::RungeKuttaBase::Solve] solution with constraints not yet implemented.");
                        // double timeDependentConstraintFactor = this->CalculateTimeDependentConstraintFactor(curTime);
                        // mStructure->ConstraintSetRHS(mTimeDependentConstraint,timeDependentConstraintFactor);
                        // mStructure->ConstraintGetRHSAfterGaussElimination(bRHS);
                    }
                    // calculate external force
                    extLoad = CalculateCurrentExternalLoad(curTime);
                }

                dof_dt0_tmp.K = mStructure->NodeCalculateDependentDofValues(dof_dt0_tmp.J);


                // velocities update
                // dof_dt1_tmp.K = mStructure->NodeCalculateDependentDofValues(dof_dt1_tmp.J);


                mStructure->NodeMergeDofValues(0, dof_dt0_tmp);

                // merge velocities
                // mStructure->NodeMergeDofValues(1,dof_dt1_tmp);

                mStructure->ElementTotalUpdateTmpStaticData();

                // calculate internal force (with update of history variables = true)
                intForce = mStructure->BuildGlobalInternalGradient();

                // ===> constraints
                intForce.ApplyCMatrix(mStructure->GetConstraintMatrix());

                // update derivatives (ydot or k1,k2,k3,k4) for Runge Kutta
                d_dof_dt0_tmp[countStage] = dof_dt1_tmp * mTimeStep;
// std::cout << "d_disp_j_tmp " << d_disp_j_tmp[countStage](0) << std::endl;

// ===> constraints
// d_dof_dt1_tmp[countStage]  = hessian2*(extLoad-intForce)*mTimeStep;

#if defined(HAVE_PARDISO) && defined(_OPENMP)
                NuTo::FullVector<double, Eigen::Dynamic> resultForSolver;
                mySolver.Solution(matrixForSolver, ((extLoad.J - intForce.J) * mTimeStep).Export(), resultForSolver);
                d_dof_dt1_tmp[countStage].J =
                        BlockFullVector<double>(-resultForSolver, this->mStructure->GetDofStatus());
#else
                NuTo::FullVector<double, Eigen::Dynamic> resultForSolver;
                mySolver.Solution(((extLoad.J - intForce.J) * mTimeStep).Export(), resultForSolver);
                d_dof_dt1_tmp[countStage].J = BlockFullVector<double>(
                        -resultForSolver, this->mStructure->GetDofStatus()); // ?? true as parameter
#endif
                // d_dof_dt1_tmp[countStage].J = -1*mStructure->SolveBlockSystem(hessian2.JJ, (extLoad.J -
                // intForce.J)*mTimeStep);

                // std::cout << "d_vel_j_tmp " << d_vel_j_tmp[countStage](0) << std::endl;
                // std::cout << "norm of acc " << (d_vel_j_tmp).norm() << std::endl;

                dof_dt0_new.J += d_dof_dt0_tmp[countStage].J * (GetStageWeights(countStage));
                // std::cout << "disp_j_new " << disp_j_new(0) << std::endl;
                dof_dt1_new.J += d_dof_dt1_tmp[countStage].J * (GetStageWeights(countStage));
                // std::cout << "vel_j_new " << vel_j_new(0) << std::endl;
            }

            mTime = prevTime + mTimeStep;
            curTime = prevCurTime + mTimeStep;

            // std::cout << "final disp_j_new " << disp_j_new(0) << std::endl;
            dof_dt0_new.K = mStructure->NodeCalculateDependentDofValues(dof_dt0_new.J);

            // velocities
            // dof_dt1_new.K = mStructure->NodeCalculateDependentDofValues(dof_dt1_new.J);

            mStructure->NodeMergeDofValues(0, dof_dt0_new);


            //            mStructure->NodeMergeDofValues(1, dof_dt1_new);

            mStructure->ElementTotalUpdateTmpStaticData();
            mStructure->ElementTotalUpdateStaticData();
            // std::cout << "delta disp between time steps" <<  (disp_j-disp_j_new).norm() << std::endl;
            dof_dt0 = dof_dt0_new;
            dof_dt1 = dof_dt1_new;

            //**********************************************
            // PostProcessing
            //**********************************************
            // outOfBalance_j is automatically zero
            // outOfBalance_j.Resize(intForce_j.GetNumRows());
            // the acceleration of the dofs k is given by the acceleration of the rhs of the constraint equation
            // this is calculated using finite differencs
            // make sure to recalculate the internal force and external force (if time factor is not 1)
            if (mTimeDependentConstraint != -1)
            {
                throw MechanicsException(
                        "[NuTo::RungeKuttaBase::Solve] solution with constraints not yet implemented.");
            }

            // acc_k = (bRHSprev-bRHShalf*2+bRHSend)*(4./(timeStep*timeStep))
            // outOfBalance_k = intForce_k - extForce_k + massMatrix_k.asDiagonal()*acc_k;

            // postprocess data for plotting
            this->PostProcess(extLoad - intForce);
        }
#if defined(HAVE_PARDISO) && defined(_OPENMP)
        mySolver.CleanUp(matrixForSolver);
#endif
    }
    catch (MechanicsException& e)
    {
        e.AddMessage("[NuTo::RungeKuttaBase::Solve] performing Newton-Raphson iteration.");
        throw e;
    }

    return NuTo::eError::SUCCESSFUL;
}


#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::RungeKuttaBase)
#endif // SWIG
#endif // ENABLE_SERIALIZATION
