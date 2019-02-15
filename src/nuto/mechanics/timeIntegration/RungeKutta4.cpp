// $Id: RungeKutta4.cpp 575 2011-09-20 18:05:35Z unger3 $

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

#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/nodes/NodeDof.h"
#include "nuto/mechanics/nodes/NodeEnum.h"
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/timeIntegration/RungeKutta4.h"
#include "nuto/mechanics/timeIntegration/TimeIntegrationEnum.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/structures/StructureBaseEnum.h"
#include "nuto/mechanics/dofSubMatrixStorage/BlockScalar.h"

#include "nuto/math/SparseMatrixCSRVector2.h"
#include "nuto/base/ErrorEnum.h"
#include "nuto/mechanics/structures/StructureOutputBlockMatrix.h"
#include "nuto/base/Timer.h"
#include <eigen3/Eigen/Dense>

//! @brief constructor
//! @param mDimension number of nodes
NuTo::RungeKutta4::RungeKutta4(StructureBase* rStructure)
    : RungeKuttaBase(rStructure)
{
}


//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::RungeKutta4::Info() const
{
    TimeIntegrationBase::Info();
}

//! @brief calculate the critical time step for explicit routines
//! for implicit routines, this will simply return zero (cmp HasCriticalTimeStep())
//! this is the critical time step from velocity verlet, the real one is certainly larger
double NuTo::RungeKutta4::CalculateCriticalTimeStep() const
{
    double maxGlobalEigenValue = mStructure->ElementTotalCalculateLargestElementEigenvalue();
    return 2.8 / std::sqrt(maxGlobalEigenValue);
}

//! @brief ... return delta time factor of intermediate stages (c in Butcher tableau, but only the delta to the previous
//! step)
// so essentially it's c_n-c_(n-1)
double NuTo::RungeKutta4::GetStageTimeFactor(int rStage) const
{
    assert(rStage < 4);
    double s;
    switch (rStage)
    {
    case 0:
        s = 0.;
        break;
    case 1:
        s = 0.5;
        break;
    case 2:
        s = 0.5;
        break;
    case 3:
        s = 1.0;
        break;
    default:
        throw MechanicsException("[NuTo::RungeKutta4::GetStageTimeFactor] rStage>3 not implemented.");
    }
    return s;
}

//! @brief ... return delta time factor of intermediate stages (c in Butcher tableau, but only the delta to the previous
//! step)
// so essentially it's c_n-c_(n-1)
bool NuTo::RungeKutta4::HasTimeChanged(int rStage) const
{
    assert(rStage < 4);
    bool s;
    switch (rStage)
    {
    case 0:
        s = false; // same as last step from the last iteration
        break;
    case 1:
        s = true;
        break;
    case 2:
        s = false;
        break;
    case 3:
        s = true;
        break;
    default:
        throw MechanicsException("[NuTo::RungeKutta4::HasTimeChanged] rStage>3 not implemented.");
    }
    return s;
}


//! @brief ... return scaling for the intermediate stage for y (a in Butcher tableau)
void NuTo::RungeKutta4::GetStageDerivativeFactor(std::vector<double>& rWeight, int rStage) const
{
    assert(rStage < 4);
    assert(rWeight.size() == 3);
    switch (rStage)
    {
    case 0:
        break;
    case 1:
        rWeight[0] = 0.5;
        break;
    case 2:
        rWeight[0] = 0.0;
        rWeight[1] = 0.5;
        break;
    case 3:
        rWeight[0] = 0.0;
        rWeight[1] = 0.0;
        rWeight[2] = 1.0;
        break;
    default:
        throw MechanicsException("[NuTo::RungeKutta4::GetStageDerivativeFactor] rStage>3 not implemented.");
    }
}

//! @brief ... return weights for the intermediate stage for y (b in Butcher tableau)
double NuTo::RungeKutta4::GetStageWeights(int rStage) const
{
    assert(rStage < 4);
    double s;
    switch (rStage)
    {
    case 0:
        s = 1. / 6.;
        break;
    case 1:
        s = 1. / 3.;
        break;
    case 2:
        s = 1. / 3.;
        break;
    case 3:
        s = 1. / 6.;
        break;
    default:
        throw MechanicsException("[NuTo::RungeKutta4::GetStageWeights] rStage>3 not implemented.");
    }
    return s;
}

//! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the
//! file      in case of restoring from a file with the wrong object type, the file id is printed
//! @return    class name
std::string NuTo::RungeKutta4::GetTypeId() const
{
    return std::string("RungeKutta4");
}

void NuTo::RungeKutta4::f(NuTo::StructureBase* mStructure, NuTo::SparseDirectSolverMUMPS& mySolver,
                          const NuTo::StructureOutputBlockVector& extLoad,
                          const NuTo::StructureOutputBlockVector& dof_dt0,
                          const NuTo::StructureOutputBlockVector& dof_dt1, double factor,
                          const NuTo::StructureOutputBlockVector& acceleration0,
                          const NuTo::StructureOutputBlockVector& velocity0,
                          NuTo::StructureOutputBlockVector& rAcceleration, NuTo::StructureOutputBlockVector& rVelocity)
{
    NuTo::StructureOutputBlockVector dof_dt0_temp(mStructure->GetDofStatus(), true);
    NuTo::StructureOutputBlockVector dof_dt1_temp(mStructure->GetDofStatus(), true);

    NuTo::StructureOutputBlockVector residual(mStructure->GetDofStatus(), true);
    NuTo::BlockFullVector<double> residual_mod(mStructure->GetDofStatus());
    const auto& cmat = mStructure->GetConstraintMatrix();

    //    std::cout << "    f before Inf norm vel: " << velocity0.J.CalculateInfNorm() << std::endl;
    //    std::cout << "    f before Inf norm acc: " << acceleration0.J.CalculateInfNorm() << std::endl;

    dof_dt0_temp.J = dof_dt0.J + factor * mTimeStep * velocity0.J; // u
    dof_dt1_temp.J = dof_dt1.J + factor * mTimeStep * acceleration0.J; // v

    dof_dt0_temp.K = mStructure->NodeCalculateDependentDofValues(dof_dt0_temp.J);
    mStructure->NodeMergeDofValues(0, dof_dt0_temp);

    dof_dt1_temp.K = mStructure->NodeCalculateDependentDofValues(dof_dt1_temp.J);
    mStructure->NodeMergeDofValues(1, dof_dt1_temp);

    //    NuTo::StructureOutputBlockVector displacements = mStructure->NodeExtractDofValues(0);
    //    std::cout << "    f before Inf norm displacements: " << displacements.J.CalculateInfNorm() << std::endl;

    NuTo::StructureOutputBlockVector intForce = mStructure->BuildGlobalInternalGradient();
    //    std::cout << "    f before Inf norm intforce: " << intForce.J.CalculateInfNorm() << std::endl;

    residual = extLoad - intForce;
    residual.ApplyCMatrix(residual_mod, cmat);

    NuTo::FullVector<double, Eigen::Dynamic> resultForSolver;
    mySolver.Solution(residual_mod.Export(), resultForSolver);
    rAcceleration.J = NuTo::BlockFullVector<double>(resultForSolver, mStructure->GetDofStatus());
    rAcceleration.K = mStructure->NodeCalculateDependentDofValues(rAcceleration.J);

    rVelocity = mStructure->NodeExtractDofValues(1);

    //    std::cout << "    f after Inf norm vel: " << rVelocity.J.CalculateInfNorm() << std::endl;
    //    std::cout << "    f after Inf norm acc: " << rAcceleration.J.CalculateInfNorm() << std::endl;
    //    std::cout << "============================================" << std::endl;

    //    mWriteForced = true;
    //    this->PostProcess(extLoad);
    //    mWriteForced = false;
}


void NuTo::RungeKutta4::f_for1Dcontact(
        NuTo::StructureBase* mStructure, const std::function<double(double)>& penaltyLaw, int dof1Dcontact,
        NuTo::SparseDirectSolverMUMPS& mySolver, const NuTo::StructureOutputBlockVector& extLoad,
        const NuTo::StructureOutputBlockVector& dof_dt0, const NuTo::StructureOutputBlockVector& dof_dt1, double factor,
        const NuTo::StructureOutputBlockVector& acceleration0, const NuTo::StructureOutputBlockVector& velocity0,
        NuTo::StructureOutputBlockVector& rAcceleration, NuTo::StructureOutputBlockVector& rVelocity)
{
    NuTo::StructureOutputBlockVector dof_dt0_temp(mStructure->GetDofStatus(), true);
    NuTo::StructureOutputBlockVector dof_dt1_temp(mStructure->GetDofStatus(), true);

    NuTo::StructureOutputBlockVector residual(mStructure->GetDofStatus(), true);
    NuTo::BlockFullVector<double> residual_mod(mStructure->GetDofStatus());
    const auto& cmat = mStructure->GetConstraintMatrix();

    //    std::cout << "    f before Inf norm vel: " << velocity0.J.CalculateInfNorm() << std::endl;
    //    std::cout << "    f before Inf norm acc: " << acceleration0.J.CalculateInfNorm() << std::endl;

    dof_dt0_temp.J = dof_dt0.J + factor * mTimeStep * velocity0.J; // u
    dof_dt1_temp.J = dof_dt1.J + factor * mTimeStep * acceleration0.J; // v

    dof_dt0_temp.K = mStructure->NodeCalculateDependentDofValues(dof_dt0_temp.J);
    mStructure->NodeMergeDofValues(0, dof_dt0_temp);

    dof_dt1_temp.K = mStructure->NodeCalculateDependentDofValues(dof_dt1_temp.J);
    mStructure->NodeMergeDofValues(1, dof_dt1_temp);

    //    NuTo::StructureOutputBlockVector displacements = mStructure->NodeExtractDofValues(0);
    //    std::cout << "    f before Inf norm displacements: " << displacements.J.CalculateInfNorm() << std::endl;

    NuTo::StructureOutputBlockVector intForce = mStructure->BuildGlobalInternalGradient();
    //    std::cout << "    f before Inf norm intforce: " << intForce.J.CalculateInfNorm() << std::endl;

    double dispIncontact = dof_dt0_temp.J[NuTo::Node::eDof::DISPLACEMENTS](dof1Dcontact);
    if (dispIncontact < 0)
        intForce.J[NuTo::Node::eDof::DISPLACEMENTS](dof1Dcontact) += penaltyLaw(dispIncontact);

    residual = extLoad - intForce;
    residual.ApplyCMatrix(residual_mod, cmat);

    NuTo::FullVector<double, Eigen::Dynamic> resultForSolver;
    mySolver.Solution(residual_mod.Export(), resultForSolver);
    rAcceleration.J = NuTo::BlockFullVector<double>(resultForSolver, mStructure->GetDofStatus());
    rAcceleration.K = mStructure->NodeCalculateDependentDofValues(rAcceleration.J);

    rVelocity = mStructure->NodeExtractDofValues(1);

    //    std::cout << "    f after Inf norm vel: " << rVelocity.J.CalculateInfNorm() << std::endl;
    //    std::cout << "    f after Inf norm acc: " << rAcceleration.J.CalculateInfNorm() << std::endl;
    //    std::cout << "============================================" << std::endl;

    //    mWriteForced = true;
    //    this->PostProcess(extLoad);
    //    mWriteForced = false;
}


//! @brief perform the time integration
//! @param rStructure ... structure
//! @param rTimeDelta ... length of the simulation
NuTo::eError NuTo::RungeKutta4::SolveRK4(double rTimeDelta, int rLoadCase)
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

        CalculateStaticAndTimeDependentExternalLoad();

        NuTo::StructureOutputBlockMatrix hessian2 = mStructure->BuildGlobalHessian2Lumped();

        hessian2.ApplyCMatrix(mStructure->GetConstraintMatrix());

        // allocate solver
        NuTo::SparseDirectSolverMUMPS mySolver;
        const BlockSparseMatrix& rMatrix = hessian2.JJ;
        NuTo::SparseMatrixCSRGeneral<double> matrixForSolver = rMatrix.ExportToCSRGeneral();
        matrixForSolver.SetOneBasedIndexing();
        mySolver.Factorization(matrixForSolver);

        double curTime = 0;

        NuTo::StructureOutputBlockVector dof_dt0 = mStructure->NodeExtractDofValues(0);
        NuTo::StructureOutputBlockVector dof_dt1 = mStructure->NodeExtractDofValues(1);

        std::vector<StructureOutputBlockVector> kAcc(this->GetNumStages() + 1, mStructure->GetDofStatus());
        std::vector<StructureOutputBlockVector> kVel(this->GetNumStages() + 1, mStructure->GetDofStatus());

        kAcc[0] = mStructure->NodeExtractDofValues(0); // displacements for all dofs
        kVel[0] = mStructure->NodeExtractDofValues(0); // displacements for all dofs

        kAcc[0].SetZero();
        kVel[0].SetZero();
        // the external load is not deformation dependent
        StructureOutputBlockVector extLoad(mStructure->GetDofStatus(), true);

        mWriteForced = true;
        this->PostProcess(extLoad);
        mWriteForced = false;

        // hessian2.CwiseInvert();
        // std::cout << "Inverse mass: \n" << hessian2.JJ.ExportToFullMatrix() << std::endl;

        while (curTime < rTimeDelta)
        {
            std::cout << "==>curTime " << curTime << " (" << curTime / rTimeDelta
                      << ") max Disp = " << dof_dt0.J[Node::eDof::DISPLACEMENTS].cwiseAbs().maxCoeff() << std::endl;

            for (int i = 0; i < this->GetNumStages(); i++)
            {
                extLoad = mStructure->BuildGlobalExternalLoadVector(rLoadCase);
                extLoad.ApplyCMatrix(mStructure->GetConstraintMatrix());

                f(mStructure, mySolver, extLoad, dof_dt0, dof_dt1, this->GetStageTimeFactor(i), kAcc[i], kVel[i],
                  kAcc[i + 1], kVel[i + 1]);
            }

            // update
            for (int i = 1; i < this->GetNumStages() + 1; i++)
            {
                dof_dt0.J += mTimeStep * GetStageWeights(i - 1) * kVel[i].J;
                dof_dt1.J += mTimeStep * GetStageWeights(i - 1) * kAcc[i].J;
            }

            dof_dt0.K = mStructure->NodeCalculateDependentDofValues(dof_dt0.J);
            mStructure->NodeMergeDofValues(0, dof_dt0);

            dof_dt1.K = mStructure->NodeCalculateDependentDofValues(dof_dt1.J);
            mStructure->NodeMergeDofValues(1, dof_dt1);

            // postprocess data for plotting
            // intForce = mStructure->BuildGlobalInternalGradient();
            // this->PostProcess(extLoad - intForce);
            // !! outofbalance only needed for  GROUP_NODE_FORCE !!
            mWriteForced = true;
            this->PostProcess(extLoad);

            mTime += mTimeStep;
            curTime += mTimeStep;
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

//! @brief perform the time integration
//! @param rStructure ... structure
//! @param rTimeDelta ... length of the simulation
NuTo::eError NuTo::RungeKutta4::RK4_DoStep_1DContact(int dof1DContact, const std::function<double(double)>& penaltyLaw,
                                                     double curTime, NuTo::SparseDirectSolverMUMPS& mySolver,
                                                     std::vector<StructureOutputBlockVector>& kAcc,
                                                     std::vector<StructureOutputBlockVector>& kVel,
                                                     StructureOutputBlockVector& extLoad,
                                                     NuTo::StructureOutputBlockVector& dof_dt0,
                                                     NuTo::StructureOutputBlockVector& dof_dt1, double simulationTime)
{
    NuTo::Timer timer(__PRETTY_FUNCTION__, mStructure->GetShowTime(), mStructure->GetLogger());

    try
    {
        if (1 || dof_dt0.J[Node::eDof::DISPLACEMENTS].cwiseAbs().maxCoeff() > 10)
        {
            std::cout << "==>curTime " << curTime << " (" << curTime / simulationTime
                      << ") max Disp = " << dof_dt0.J[Node::eDof::DISPLACEMENTS].cwiseAbs().maxCoeff() << std::endl;
        }

        for (int i = 0; i < this->GetNumStages(); i++)
        {
            // no time dependent load
            // extLoad = mStructure->BuildGlobalExternalLoadVector(rLoadCase);
            // extLoad.ApplyCMatrix(mStructure->GetConstraintMatrix());

            f_for1Dcontact(mStructure, penaltyLaw, dof1DContact, mySolver, extLoad, dof_dt0, dof_dt1,
                           this->GetStageTimeFactor(i), kAcc[i], kVel[i], kAcc[i + 1], kVel[i + 1]);
        }

        // update
        for (int i = 1; i < this->GetNumStages() + 1; i++)
        {
            //            std::cout << "Nach f Inf norm: " << kVel[i].J.CalculateInfNorm() << std::endl;
            dof_dt0.J += mTimeStep * GetStageWeights(i - 1) * kVel[i].J;
            dof_dt1.J += mTimeStep * GetStageWeights(i - 1) * kAcc[i].J;
        }

        dof_dt0.K = mStructure->NodeCalculateDependentDofValues(dof_dt0.J);
        mStructure->NodeMergeDofValues(0, dof_dt0);

        dof_dt1.K = mStructure->NodeCalculateDependentDofValues(dof_dt1.J);
        mStructure->NodeMergeDofValues(1, dof_dt1);

        // postprocess data for plotting
        // intForce = mStructure->BuildGlobalInternalGradient();
        // this->PostProcess(extLoad - intForce);
        // !! outofbalance only needed for  GROUP_NODE_FORCE !!
        mWriteForced = true;
        this->PostProcess(extLoad);

        mTime += mTimeStep;
        curTime += mTimeStep;
    }
    catch (MechanicsException& e)
    {
        e.AddMessage("[NuTo::RungeKuttaBase::Solve] performing Newton-Raphson iteration.");
        throw e;
    }

    return NuTo::eError::SUCCESSFUL;
}


//! @brief perform the time integration
//! @param rStructure ... structure
//! @param rTimeDelta ... length of the simulation
NuTo::eError NuTo::RungeKutta4::RK4_DoStep(int rLoadCase, double curTime, NuTo::SparseDirectSolverMUMPS& mySolver,
                                           std::vector<StructureOutputBlockVector>& kAcc,
                                           std::vector<StructureOutputBlockVector>& kVel,
                                           StructureOutputBlockVector& extLoad,
                                           NuTo::StructureOutputBlockVector& dof_dt0,
                                           NuTo::StructureOutputBlockVector& dof_dt1, double simulationTime)
{
    NuTo::Timer timer(__PRETTY_FUNCTION__, mStructure->GetShowTime(), mStructure->GetLogger());

    try
    {
        if (1 || dof_dt0.J[Node::eDof::DISPLACEMENTS].cwiseAbs().maxCoeff() > 10)
        {
            std::cout << "==>curTime " << curTime << " (" << curTime / simulationTime
                      << ") max Disp = " << dof_dt0.J[Node::eDof::DISPLACEMENTS].cwiseAbs().maxCoeff() << std::endl;
        }

        for (int i = 0; i < this->GetNumStages(); i++)
        {
            // no time dependent load
            // extLoad = mStructure->BuildGlobalExternalLoadVector(rLoadCase);
            // extLoad.ApplyCMatrix(mStructure->GetConstraintMatrix());

            f(mStructure, mySolver, extLoad, dof_dt0, dof_dt1, this->GetStageTimeFactor(i), kAcc[i], kVel[i],
              kAcc[i + 1], kVel[i + 1]);
        }

        // update
        for (int i = 1; i < this->GetNumStages() + 1; i++)
        {
            //            std::cout << "Nach f Inf norm: " << kVel[i].J.CalculateInfNorm() << std::endl;
            dof_dt0.J += mTimeStep * GetStageWeights(i - 1) * kVel[i].J;
            dof_dt1.J += mTimeStep * GetStageWeights(i - 1) * kAcc[i].J;
        }

        dof_dt0.K = mStructure->NodeCalculateDependentDofValues(dof_dt0.J);
        mStructure->NodeMergeDofValues(0, dof_dt0);

        dof_dt1.K = mStructure->NodeCalculateDependentDofValues(dof_dt1.J);
        mStructure->NodeMergeDofValues(1, dof_dt1);

        // postprocess data for plotting
        // intForce = mStructure->BuildGlobalInternalGradient();
        // this->PostProcess(extLoad - intForce);
        // !! outofbalance only needed for  GROUP_NODE_FORCE !!
        mWriteForced = true;
        this->PostProcess(extLoad);

        mTime += mTimeStep;
        curTime += mTimeStep;
    }
    catch (MechanicsException& e)
    {
        e.AddMessage("[NuTo::RungeKuttaBase::Solve] performing Newton-Raphson iteration.");
        throw e;
    }

    return NuTo::eError::SUCCESSFUL;
}
