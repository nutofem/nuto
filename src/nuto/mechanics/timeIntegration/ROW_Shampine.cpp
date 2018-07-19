#ifdef _OPENMP
#include <omp.h>
#endif

#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/nodes/NodeDof.h"
#include "nuto/mechanics/nodes/NodeEnum.h"
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/timeIntegration/ROW_Shampine.h"
#include "nuto/mechanics/timeIntegration/TimeIntegrationEnum.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/structures/StructureBaseEnum.h"

#include "nuto/math/SparseMatrixCSRVector2.h"
#include "nuto/base/ErrorEnum.h"
#include "nuto/mechanics/structures/StructureOutputBlockMatrix.h"
#include "nuto/base/Timer.h"
#include <eigen3/Eigen/Dense>

//! @brief constructor
//! @param mDimension number of nodes
NuTo::ROW_Shampine::ROW_Shampine(StructureBase* rStructure)
    : RungeKuttaBase(rStructure)
{
}


//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::ROW_Shampine::Info() const
{
    TimeIntegrationBase::Info();
}

//! @brief calculate the critical time step for explicit routines
//! for implicit routines, this will simply return zero (cmp HasCriticalTimeStep())
//! this is the critical time step from velocity verlet, the real one is certainly larger
double NuTo::ROW_Shampine::CalculateCriticalTimeStep() const
{
    double maxGlobalEigenValue = mStructure->ElementTotalCalculateLargestElementEigenvalue();
    return 2.8 / std::sqrt(maxGlobalEigenValue);
}

//! @brief ... return delta time factor of intermediate stages (c in Butcher tableau, but only the delta to the previous
//! step)
// so essentially it's c_n-c_(n-1)
double NuTo::ROW_Shampine::GetStageTimeFactor(int rStage) const
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
        throw MechanicsException("[NuTo::ROW_Shampine::GetStageTimeFactor] rStage>3 not implemented.");
    }
    return s;
}

//! @brief ... return delta time factor of intermediate stages (c in Butcher tableau, but only the delta to the previous
//! step)
// so essentially it's c_n-c_(n-1)
bool NuTo::ROW_Shampine::HasTimeChanged(int rStage) const
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
        throw MechanicsException("[NuTo::ROW_Shampine::HasTimeChanged] rStage>3 not implemented.");
    }
    return s;
}


//! @brief ... return scaling for the intermediate stage for y (a in Butcher tableau)
void NuTo::ROW_Shampine::GetStageDerivativeFactor(std::vector<double>& rWeight, int rStage) const
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
        throw MechanicsException("[NuTo::ROW_Shampine::GetStageDerivativeFactor] rStage>3 not implemented.");
    }
}

//! @brief ... return weights for the intermediate stage for y (b in Butcher tableau)
double NuTo::ROW_Shampine::GetStageWeights(int rStage) const
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
        throw MechanicsException("[NuTo::ROW_Shampine::GetStageWeights] rStage>3 not implemented.");
    }
    return s;
}

//! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the
//! file      in case of restoring from a file with the wrong object type, the file id is printed
//! @return    class name
std::string NuTo::ROW_Shampine::GetTypeId() const
{
    return std::string("ROW_Shampine");
}

//! @brief perform the time integration
//! @param rStructure ... structure
//! @param rTimeDelta ... length of the simulation
NuTo::eError NuTo::ROW_Shampine::ROW_Shampine_DoStep(double rTimeDelta, double curTime,
                                                     NuTo::SparseMatrixCSRGeneral<double>& M_JJ,
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

        NuTo::StructureOutputBlockVector dof_dt0_temp(mStructure->GetDofStatus(), true);
        NuTo::StructureOutputBlockVector dof_dt1_temp(mStructure->GetDofStatus(), true);

        NuTo::StructureOutputBlockVector intForce(mStructure->GetDofStatus(), true);
        NuTo::StructureOutputBlockVector residual(mStructure->GetDofStatus(), true);
        NuTo::FullVector<double> rhs;
        NuTo::BlockFullVector<double> residual_mod(mStructure->GetDofStatus());
        NuTo::BlockFullVector<double> r1_block(mStructure->GetDofStatus());
        NuTo::FullVector<double> r1;
        NuTo::FullVector<double> r2;

        const auto& cmat = mStructure->GetConstraintMatrix();

        NuTo::StructureOutputBlockMatrix K = mStructure->BuildGlobalHessian0();
        K.ApplyCMatrix(cmat);

        double factor_acc = (2. / rTimeDelta);
        double factor_vel = (4. / (rTimeDelta * rTimeDelta));

        NuTo::SparseMatrixCSRGeneral<double> matrixForSolver = (M_JJ * factor_vel) + K.JJ.ExportToCSRGeneral();

        NuTo::SparseDirectSolverMUMPS mySolver;
        matrixForSolver.SetOneBasedIndexing();
        mySolver.Factorization(matrixForSolver);
        NuTo::FullVector<double, Eigen::Dynamic> solution;


        //%%%%%%%%%%%%%%%%%%%%%%%%%% Stage I %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
        dof_dt0_temp.J = dof_dt0.J;
        dof_dt1_temp.J = dof_dt1.J;

        dof_dt0_temp.K = mStructure->NodeCalculateDependentDofValues(dof_dt0_temp.J);
        mStructure->NodeMergeDofValues(0, dof_dt0_temp);

        dof_dt1_temp.K = mStructure->NodeCalculateDependentDofValues(dof_dt1_temp.J);
        mStructure->NodeMergeDofValues(1, dof_dt1_temp);

        intForce = mStructure->BuildGlobalInternalGradient();

        residual = extLoad - intForce;
        residual.ApplyCMatrix(residual_mod, cmat);

        rhs = factor_acc * (residual_mod.Export()) + factor_vel * (M_JJ * (dof_dt1_temp.J.Export()));

        mySolver.Solution(rhs, solution);

        kVel[0].J = NuTo::BlockFullVector<double>(solution, mStructure->GetDofStatus());
        kVel[0].K = mStructure->NodeCalculateDependentDofValues(kVel[0].J);

        kAcc[0].J = factor_acc * (kVel[0].J - dof_dt1_temp.J);
        kAcc[0].K = mStructure->NodeCalculateDependentDofValues(kAcc[0].J);

        //%%%%%%%%%%%%%%%%%%%%%%%%%% Stage II %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
        dof_dt0_temp.J = dof_dt0.J + rTimeDelta * kVel[0].J;
        dof_dt1_temp.J = dof_dt1.J + rTimeDelta * kAcc[0].J;

        dof_dt0_temp.K = mStructure->NodeCalculateDependentDofValues(dof_dt0_temp.J);
        mStructure->NodeMergeDofValues(0, dof_dt0_temp);

        dof_dt1_temp.K = mStructure->NodeCalculateDependentDofValues(dof_dt1_temp.J);
        mStructure->NodeMergeDofValues(1, dof_dt1_temp);

        intForce = mStructure->BuildGlobalInternalGradient();

        residual = extLoad - intForce;
        residual.ApplyCMatrix(residual_mod, cmat);

        r1_block = (dof_dt1_temp.J - 4 * kVel[0].J);
        r1 = r1_block.Export();
        r2 = (residual_mod - 4 * kAcc[0].J).Export();

        rhs = factor_acc * r2 + factor_vel * (M_JJ * r1);

        mySolver.Solution(rhs, solution);

        kVel[1].J = NuTo::BlockFullVector<double>(solution, mStructure->GetDofStatus());
        kVel[1].K = mStructure->NodeCalculateDependentDofValues(kVel[1].J);

        kAcc[1].J = factor_acc * (kVel[1].J - r1_block);
        kAcc[1].K = mStructure->NodeCalculateDependentDofValues(kAcc[1].J);

        //%%%%%%%%%%%%%%%%%%%%%%%%%% Stage III %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
        dof_dt0_temp.J = dof_dt0.J + (24. / 25.) * rTimeDelta * kVel[0].J + (3. / 25.) * rTimeDelta * kVel[1].J;
        dof_dt1_temp.J = dof_dt1.J + (24. / 25.) * rTimeDelta * kAcc[0].J + (3. / 25.) * rTimeDelta * kAcc[1].J;

        dof_dt0_temp.K = mStructure->NodeCalculateDependentDofValues(dof_dt0_temp.J);
        mStructure->NodeMergeDofValues(0, dof_dt0_temp);

        dof_dt1_temp.K = mStructure->NodeCalculateDependentDofValues(dof_dt1_temp.J);
        mStructure->NodeMergeDofValues(1, dof_dt1_temp);

        intForce = mStructure->BuildGlobalInternalGradient();

        residual = extLoad - intForce;
        residual.ApplyCMatrix(residual_mod, cmat);

        r1_block = dof_dt1_temp.J + (186. / 25.) * kVel[0].J + (6. / 5.) * kVel[1].J;
        r1 = r1_block.Export();
        r2 = (residual_mod + (186. / 25.) * kAcc[0].J + (6. / 5.) * kAcc[1].J).Export();

        rhs = factor_acc * r2 + factor_vel * (M_JJ * r1);

        mySolver.Solution(rhs, solution);

        kVel[2].J = NuTo::BlockFullVector<double>(solution, mStructure->GetDofStatus());
        kVel[2].K = mStructure->NodeCalculateDependentDofValues(kVel[2].J);

        kAcc[2].J = factor_acc * (kVel[2].J - r1_block);
        kAcc[2].K = mStructure->NodeCalculateDependentDofValues(kAcc[2].J);

        //%%%%%%%%%%%%%%%%%%%%%%%%%% Stage IV %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

        r1_block = dof_dt1_temp.J - (56. / 125.) * kVel[0].J - (27. / 125.) * kVel[1].J - (1. / 5.) * kVel[2].J;
        r1 = r1_block.Export();
        r2 = (residual_mod - (56. / 125.) * kAcc[0].J - (27. / 125.) * kAcc[1].J - (1. / 5.) * kAcc[2].J).Export();

        rhs = factor_acc * r2 + factor_vel * (M_JJ * r1);

        mySolver.Solution(rhs, solution);

        kVel[3].J = NuTo::BlockFullVector<double>(solution, mStructure->GetDofStatus());
        kVel[3].K = mStructure->NodeCalculateDependentDofValues(kVel[3].J);

        kAcc[3].J = factor_acc * (kVel[3].J - r1_block);
        kAcc[3].K = mStructure->NodeCalculateDependentDofValues(kAcc[3].J);

        //%%%%%%%%%%%%%%%%%%%%%%%%%% update (finale) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
        dof_dt0.J += rTimeDelta * ((19. / 18.) * kVel[0].J + (1. / 4.) * kVel[1].J + (25. / 216.) * kVel[2].J +
                                   (125. / 216.) * kVel[3].J);
        dof_dt1.J += rTimeDelta * ((19. / 18.) * kAcc[0].J + (1. / 4.) * kAcc[1].J + (25. / 216.) * kAcc[2].J +
                                   (125. / 216.) * kAcc[3].J);

        dof_dt0.K = mStructure->NodeCalculateDependentDofValues(dof_dt0.J);
        mStructure->NodeMergeDofValues(0, dof_dt0);

        dof_dt1.K = mStructure->NodeCalculateDependentDofValues(dof_dt1.J);
        mStructure->NodeMergeDofValues(1, dof_dt1);

        // postprocess data for plotting
        // intForce = mStructure->BuildGlobalInternalGradient();
        // this->PostProcess(extLoad - intForce);
        // !! outofbalance only needed for  GROUP_NODE_FORCE !!
        this->PostProcess(extLoad);

        mTime += rTimeDelta;
        curTime += rTimeDelta;
    }
    catch (MechanicsException& e)
    {
        e.AddMessage("[NuTo::RungeKuttaBase::Solve] performing Newton-Raphson iteration.");
        throw e;
    }

    return NuTo::eError::SUCCESSFUL;
}
