#include "BoostUnitTest.h"

#include "mechanics/dofs/DofNumbering.h"
#include "mechanics/cell/SimpleAssembler.h"

#include "mechanics/integrands/DynamicMomentumBalance.h"
#include "mechanics/interpolation/InterpolationTrussLobatto.h"
#include "mechanics/integrationtypes/IntegrationTypeTensorProduct.h"
#include "mechanics/constraints/Constraints.h"
#include "mechanics/mesh/UnitMeshFem.h"
#include "mechanics/mesh/MeshFemDofConvert.h"
#include "mechanics/tools/CellStorage.h"

#include "mechanics/constitutive/LinearElastic.h"

#include <Eigen/Eigenvalues>

using namespace NuTo;

struct Steel
{
    static constexpr double E = 200.0;
    static constexpr double nu = 0.3;
    static constexpr double rho = 8000.;
};

class VibratingTruss
{
public:
    VibratingTruss(int numElements, int order)
        : mMesh(UnitMeshFem::CreateLines(numElements))
        , mOrder(order)
        , mDof("Dispacement", 1)
        , mLaw(Steel::E, Steel::nu)
        , mPDE(mDof, mLaw, Steel::rho)
        , mIntegrationType(mOrder + 1, eIntegrationMethod::GAUSS)
    {
        AddDofInterpolation(&mMesh, mDof, mMesh.CreateInterpolation(InterpolationTrussLobatto(order)));

        mCellGroup = mCells.AddCells(mMesh.ElementsTotal(), mIntegrationType);
    }

    Eigen::VectorXd Solve()
    {
        DofInfo dofInfo = DofNumbering::Build(mMesh.NodesTotal(mDof), mDof, mConstraints);
        Eigen::SparseMatrix<double> cmat = mConstraints.BuildConstraintMatrix(mDof, dofInfo.numIndependentDofs[mDof]);
        SimpleAssembler asmbl = SimpleAssembler(dofInfo);

        auto stiffnessF = [&](const auto& cipd) { return mPDE.Hessian0(cipd, 0.0); };
        auto massF = [&](const auto& cipd) { return mPDE.Hessian2(cipd); };

        Eigen::SparseMatrix<double> K = asmbl.BuildMatrix(mCellGroup, {mDof}, stiffnessF).JJ(mDof, mDof);
        Eigen::SparseMatrix<double> M = asmbl.BuildMatrix(mCellGroup, {mDof}, massF).JJ(mDof, mDof);

        return Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd>(Eigen::MatrixXd(K), Eigen::MatrixXd(M))
                .eigenvalues();
    }

    Eigen::VectorXd Expected(int n)
    {
        Eigen::VectorXd result(n);
        for (int i = 0; i < n; i++)
        {
            result[i] = (i * M_PI) * (i * M_PI) * Steel::E / Steel::rho;
        }
        return result;
    }

private:
    MeshFem mMesh;
    int mOrder;

    DofType mDof;
    Laws::LinearElastic<1> mLaw;
    Integrands::DynamicMomentumBalance<1> mPDE;
    Constraint::Constraints mConstraints;

    IntegrationTypeTensorProduct<1> mIntegrationType;
    CellStorage mCells;
    Group<CellInterface> mCellGroup;
};

//!@brief Computes the first squared eigenfrequencies of a vibrating truss with free ends
//! and compares with analytical solution. Since the true solution and the approximated one
//! will differ on smaller scales / higher frequencies only low frequencies are compared.
BOOST_AUTO_TEST_CASE(VibratingTruss1DFirstEigenfrequencies)
{
    VibratingTruss example1(10, 8);
    Eigen::VectorXd computed = example1.Solve();
    Eigen::VectorXd expected = example1.Expected(computed.size());
    for (int i = 1; i < 20; i++)
    {
        BOOST_CHECK_CLOSE(computed[i], expected[i], 1.e-3);
    }
}
