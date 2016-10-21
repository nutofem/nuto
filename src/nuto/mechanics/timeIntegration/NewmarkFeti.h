#include "nuto/mechanics/timeIntegration/TimeIntegrationBase.h"
#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"
#include "nuto/base/ErrorEnum.h"
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/base/Timer.h"

#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/structures/unstructured/StructureFETI.h"
#include "nuto/mechanics/structures/StructureBaseEnum.h"
#include "nuto/mechanics/structures/StructureOutputBlockMatrix.h"
#include "nuto/mechanics/structures/StructureOutputDummy.h"
#include <mpi.h>

#include "nuto/base/CallbackInterface.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include <eigen3/Eigen/Dense>
#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/nodes/NodeEnum.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveCalculateStaticData.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveTimeStep.h"

#include <cmath>

namespace NuTo
{
class NewmarkFeti : public NewmarkDirect
{
public:
    using VectorXd      = Eigen::VectorXd;
    using MatrixXd      = Eigen::MatrixXd;
    using SparseMatrix  = Eigen::SparseMatrix<double>;

    ///
    /// \brief NewmarkFeti
    /// \param rStructure
    ///
    NewmarkFeti(StructureBase* rStructure) : NewmarkDirect (rStructure)
    {

    }

    ///
    /// \brief GetTypeId
    /// \return
    ///
    std::string GetTypeId() const override
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented!");
    }

    ///
    /// \brief HasCriticalTimeStep
    /// \return
    ///
    bool HasCriticalTimeStep() const override
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented!");
    }

    ///
    /// \brief CalculateCriticalTimeStep
    /// \return
    ///
    double CalculateCriticalTimeStep() const override
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented!");
    }


    ///
    /// \brief GatherInterfaceRigidBodyModes
    /// \param interfaceRigidBodyModes
    /// \param numRigidBodyModesGlobal
    /// \return
    ///
    MatrixXd GatherInterfaceRigidBodyModes(const Eigen::MatrixXd& interfaceRigidBodyModes, const int numRigidBodyModesGlobal);

    ///
    /// \brief GatherRigidBodyForceVector
    /// \param rigidBodyForceVectorLocal
    /// \param numRigidBodyModesGlobal
    /// \return
    ///
    VectorXd GatherRigidBodyForceVector(const Eigen::VectorXd& rigidBodyForceVectorLocal, const int numRigidBodyModesGlobal);

    ///
    /// \brief MpiGatherRecvCountAndDispls
    /// \param recvCount
    /// \param displs
    /// \param numValues
    ///
    void MpiGatherRecvCountAndDispls(std::vector<int>& recvCount, std::vector<int>& displs, const int numValues);

    ///
    /// \brief CalculateNormResidual
    /// \param residual_mod
    /// \param activeDofSet
    /// \return
    ///
    BlockScalar CalculateNormResidual(BlockFullVector<double>& residual_mod, const std::set<Node::eDof>& activeDofSet);

    //! @brief Projected stabilized Bi-conjugate gradient method (BiCGStab)
    //!
    //! Iterative method to solve the interface problem.
    //! @param projection: Projection matrix that guarantees that the solution satisfies the constraint at every iteration
    int BiCgStab(const MatrixXd& projection, VectorXd& x, const VectorXd& rhs);

    //! @brief Conjugate projected gradient method (CG)
    //!
    //! Iterative method to solve the interface problem.
    //! @param projection: Projection matrix that guarantees that the solution satisfies the constraint at every iteration
    int CPG(const MatrixXd& projection, VectorXd& x, const VectorXd& rhs);

    //! @brief Solves the global and local problem
    //!
    //! Calculates the rigid body modes of the structure.
    //! Calculates the initial guess for the projected CG/BiCGStab method
    //! Solves for the Lagrange multipliers at the subdomain interfaces.
    //! Calculates the increment of the free degrees of freedom
    //!
    StructureOutputBlockVector FetiSolve(BlockFullVector<double> residual_mod, const std::set<Node::eDof>& activeDofSet, VectorXd& deltaLambda);

    //! @brief perform the time integration
    //! @param rTimeDelta ... length of the simulation
    NuTo::eError Solve(double rTimeDelta) override;

private:
    Eigen::SparseQR<Eigen::SparseMatrix<double>,Eigen::COLAMDOrdering<int>> mSolver;
    SparseMatrix mLocalPreconditioner;
    SparseMatrix mTangentStiffnessMatrix;
    const double    mCpgTolerance     = 1.0e-6;
    const int       mCpgMaxIterations = 1000;
};
}// namespace NuTo
