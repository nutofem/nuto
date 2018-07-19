// $Id: RungeKuttaDormandPrince.cpp 575 2011-09-20 18:05:35Z unger3 $

#ifdef _OPENMP
#include <omp.h>
#endif

#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/timeIntegration/RungeKuttaDormandPrince.h"
#include "nuto/mechanics/timeIntegration/TimeIntegrationEnum.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/nodes/NodeEnum.h"
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/structures/StructureBaseEnum.h"

#include "nuto/math/SparseMatrixCSRVector2.h"
#include "nuto/base/ErrorEnum.h"
#include "nuto/mechanics/structures/StructureOutputBlockMatrix.h"
#include "nuto/base/Timer.h"
#include <eigen3/Eigen/Dense>


//! @brief constructor
//! @param mDimension number of nodes
NuTo::RungeKuttaDormandPrince::RungeKuttaDormandPrince(StructureBase* rStructure)
    : RungeKuttaBase(rStructure)
{
}


//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::RungeKuttaDormandPrince::Info() const
{
    TimeIntegrationBase::Info();
}

//! @brief calculate the critical time step for explicit routines (this is wrong do)
//! for implicit routines, this will simply return zero (cmp HasCriticalTimeStep())
//! this is the critical time step from velocity verlet, the real one is certainly larger
double NuTo::RungeKuttaDormandPrince::CalculateCriticalTimeStep() const
{
    double maxGlobalEigenValue = mStructure->ElementTotalCalculateLargestElementEigenvalue();
    return 2.8 / std::sqrt(maxGlobalEigenValue);
}

//! @brief ... return delta time factor of intermediate stages (c in Butcher tableau, but only the delta to the previous
//! step)
// so essentially it's c_n-c_(n-1)
double NuTo::RungeKuttaDormandPrince::GetStageTimeFactor(int rStage) const
{
    assert(rStage < 7);
    double s;
    switch (rStage)
    {
    case 0:
        s = 0.;
        break;
    case 1:
        s = 0.2;
        break;
    case 2:
        s = 0.3;
        break;
    case 3:
        s = 0.8;
        break;
    case 4:
        s = 8. / 9.;
        break;
    case 5:
        s = 1.0;
        break;
    case 6:
        s = 1.0;
        break;
    default:
        throw MechanicsException("[NuTo::RungeKuttaDormandPrince::GetStageTimeFactor] rStage<7.");
    }
    return s;
}

//! @brief ... return delta time factor of intermediate stages (c in Butcher tableau, but only the delta to the previous
//! step)
// so essentially it's c_n-c_(n-1)
bool NuTo::RungeKuttaDormandPrince::HasTimeChanged(int rStage) const
{
    assert(rStage < 7);
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
        s = true;
        break;
    case 3:
        s = true;
        break;
    case 4:
        s = true;
        break;
    case 5:
        s = true;
        break;
    case 6:
        s = false;
        break;
    default:
        throw MechanicsException("[NuTo::RungeKuttaDormandPrince::HasTimeChanged] rStage<7.");
    }
    return s;
}


//! @brief ... return scaling for the intermediate stage for y (a in Butcher tableau)
void NuTo::RungeKuttaDormandPrince::GetStageDerivativeFactor(std::vector<double>& rWeight, int rStage) const
{
    assert(rStage < 7);
    assert(rWeight.size() == 6);
    switch (rStage)
    {
    case 0:
        break;
    case 1:
        rWeight[0] = 0.2;
        break;
    case 2:
        rWeight[0] = 3. / 40.;
        rWeight[1] = 9. / 40.;
        break;
    case 3:
        rWeight[0] = 44. / 45.;
        rWeight[1] = -56. / 15.;
        rWeight[2] = 32. / 9.;
        break;
    case 4:
        rWeight[0] = 19372. / 6561.;
        rWeight[1] = -25360. / 2187.;
        rWeight[2] = 64448. / 6561.;
        rWeight[3] = -212. / 729.;
        break;
    case 5:
        rWeight[0] = 9017. / 3168.;
        rWeight[1] = -355. / 33.;
        rWeight[2] = 46732. / 5247.;
        rWeight[3] = 49. / 176.;
        rWeight[4] = -5103. / 18656.;
        break;
    case 6:
        rWeight[0] = 35. / 384.;
        rWeight[1] = 0.;
        rWeight[2] = 500. / 1113.;
        rWeight[3] = 125. / 192.;
        rWeight[4] = -2187. / 6784.;
        rWeight[5] = 11. / 84.;
        break;
    default:
        throw MechanicsException("[NuTo::RungeKuttaDormandPrince::GetStageDerivativeFactor] rStage<7.");
    }
}

//! @brief ... return weights for the intermediate stage for y (b in Butcher tableau)
double NuTo::RungeKuttaDormandPrince::GetStageWeights(int rStage) const
{
    assert(rStage < 7);
    double s;
    switch (rStage)
    {
    case 0:
        s = 35. / 385;
        break;
    case 1:
        s = 0.;
        break;
    case 2:
        s = 500. / 1113.;
        break;
    case 3:
        s = 125. / 192.;
        break;
    case 4:
        s = -2187. / 6784.;
        break;
    case 5:
        s = 11. / 84.;
        break;
    case 6:
        s = 0.;
        break;
    default:
        throw MechanicsException("[NuTo::RungeKuttaDormandPrince::GetStageWeights] rStage<7.");
    }
    return s;
}

//! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the
//! file
//!            in case of restoring from a file with the wrong object type, the file id is printed
//! @return    class name
std::string NuTo::RungeKuttaDormandPrince::GetTypeId() const
{
    return std::string("RungeKuttaDormandPrince");
}

void NuTo::RungeKuttaDormandPrince::f(NuTo::StructureBase* mStructure, NuTo::SparseDirectSolverMUMPS& mySolver,
                                      const NuTo::StructureOutputBlockVector& extLoad,
                                      const NuTo::StructureOutputBlockVector& dof_dt0,
                                      const NuTo::StructureOutputBlockVector& dof_dt1, double factor,
                                      const NuTo::StructureOutputBlockVector& acceleration0,
                                      const NuTo::StructureOutputBlockVector& velocity0,
                                      NuTo::StructureOutputBlockVector& rAcceleration,
                                      NuTo::StructureOutputBlockVector& rVelocity, int stage)
{
    NuTo::StructureOutputBlockVector dof_dt0_temp(mStructure->GetDofStatus(), true);
    NuTo::StructureOutputBlockVector dof_dt1_temp(mStructure->GetDofStatus(), true);

    NuTo::StructureOutputBlockVector residual(mStructure->GetDofStatus(), true);
    NuTo::BlockFullVector<double> residual_mod(mStructure->GetDofStatus());
    const auto& cmat = mStructure->GetConstraintMatrix();

    dof_dt0_temp.J = dof_dt0.J; // u
    dof_dt1_temp.J = dof_dt1.J; // v

    std::vector<double> factors;
    factors.resize(6);

    GetStageDerivativeFactor(factors, stage);

    for (int i = 1; i < stage; i++)
    {
        dof_dt0_temp.J += factors[i] * mTimeStep * velocity0.J; // u
        dof_dt1_temp.J += factors[i] * mTimeStep * acceleration0.J; // v
    }

    dof_dt0_temp.K = mStructure->NodeCalculateDependentDofValues(dof_dt0_temp.J);
    mStructure->NodeMergeDofValues(0, dof_dt0_temp);

    dof_dt1_temp.K = mStructure->NodeCalculateDependentDofValues(dof_dt1_temp.J);
    mStructure->NodeMergeDofValues(1, dof_dt1_temp);

    NuTo::StructureOutputBlockVector intForce = mStructure->BuildGlobalInternalGradient();

    residual = extLoad - intForce;
    residual.ApplyCMatrix(residual_mod, cmat);

    NuTo::FullVector<double, Eigen::Dynamic> resultForSolver;
    mySolver.Solution(residual_mod.Export(), resultForSolver);
    rAcceleration.J = NuTo::BlockFullVector<double>(resultForSolver, mStructure->GetDofStatus());
    rAcceleration.K = mStructure->NodeCalculateDependentDofValues(rAcceleration.J);

    rVelocity = mStructure->NodeExtractDofValues(1);
}

//! @brief perform the time integration
//! @param rStructure ... structure
//! @param rTimeDelta ... length of the simulation
NuTo::eError NuTo::RungeKuttaDormandPrince::DP_DoStep(double rTimeDelta, int rLoadCase, double curTime,
                                                      NuTo::SparseDirectSolverMUMPS& mySolver,
                                                      std::vector<StructureOutputBlockVector>& kAcc,
                                                      std::vector<StructureOutputBlockVector>& kVel,
                                                      StructureOutputBlockVector& extLoad,
                                                      NuTo::StructureOutputBlockVector& dof_dt0,
                                                      NuTo::StructureOutputBlockVector& dof_dt1)
{
    NuTo::Timer timer(__PRETTY_FUNCTION__, mStructure->GetShowTime(), mStructure->GetLogger());

    try
    {
        std::cout << "==>curTime " << curTime << " (" << curTime / rTimeDelta
                  << ") max Disp = " << dof_dt0.J[Node::eDof::DISPLACEMENTS].cwiseAbs().maxCoeff() << std::endl;

        for (int i = 0; i < this->GetNumStages(); i++)
        {
            extLoad = mStructure->BuildGlobalExternalLoadVector(rLoadCase);
            extLoad.ApplyCMatrix(mStructure->GetConstraintMatrix());

            f(mStructure, mySolver, extLoad, dof_dt0, dof_dt1, this->GetStageTimeFactor(i), kAcc[i], kVel[i],
              kAcc[i + 1], kVel[i + 1], i);
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
