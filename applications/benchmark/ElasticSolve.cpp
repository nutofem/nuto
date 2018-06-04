#include <benchmark/benchmark.h>

/*
 * Benchmarks a linear elastic linear 3d brick structure. May serve as a comparision with old nuto or other tools.
 *
 */

#include "nuto/mechanics/integrands/MomentumBalance.h"
#include "nuto/mechanics/constitutive/LinearElastic.h"

#include "nuto/mechanics/mesh/UnitMeshFem.h"
#include "nuto/mechanics/mesh/MeshFemDofConvert.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeTensorProduct.h"
#include "nuto/mechanics/tools/CellStorage.h"
#include "nuto/mechanics/tools/TimeDependentProblem.h"
#include "nuto/mechanics/tools/QuasistaticSolver.h"
#include "nuto/mechanics/constraints/ConstraintCompanion.h"

#include "nuto/mechanics/dofs/DofVectorConvertEigen.h"

using namespace NuTo;

class TestStructure
{
public:
    TestStructure(Eigen::Vector3i numElements)
        : mMesh(UnitMeshFem::CreateBricks(numElements[0], numElements[1], numElements[2]))
        , mDof("Displacement", 3)
        , mElasticLaw(30000., 0.2)
        , mMomentumBalance(mDof, mElasticLaw)
        , mIntegration(2, eIntegrationMethod::GAUSS)
        , mFunctions(&mMesh)
        , mImplicitCallBack(mFunctions, mReducedSolutionSpaceOperator)
    {
        AddDofInterpolation(&mMesh, mDof);
        Group<CellInterface> cells = mCells.AddCells(mMesh.ElementsTotal(), mIntegration);

        Constraint::Constraints constraints;
        constraints.Add(mDof, Constraint::Component(mMesh.NodesAtAxis(eDirection::X, mDof),
                                                    {eDirection::X, eDirection::Y, eDirection::Z}));
        constraints.Add(mDof, Constraint::Component(mMesh.NodesAtAxis(eDirection::X, mDof, 1), {eDirection::X},
                                                    Constraint::RhsRamp(1, 0.2)));

        mFunctions.AddGradientFunction(
                cells, TimeDependentProblem::Bind_dt(mMomentumBalance, &Integrands::MomentumBalance<3>::Gradient));
        mFunctions.AddHessian0Function(
                cells, TimeDependentProblem::Bind_dt(mMomentumBalance, &Integrands::MomentumBalance<3>::Hessian0));

        std::vector<DofType> dofTypes;
        dofTypes.push_back(mDof);
        DofContainer<int> numTotalDofs;
        DofInfo dofInfo = DofNumbering::Build(mMesh.NodesTotal(mDof), mDof, constraints);
        numTotalDofs.Insert(mDof, dofInfo.numDependentDofs[mDof] + dofInfo.numIndependentDofs[mDof]);
        mReducedSolutionSpaceOperator = ReducedSolutionSpace(dofTypes, numTotalDofs, constraints);

        DofVector<double> X = mFunctions.RenumberDofs(constraints, dofTypes, DofVector<double>());

        mSolutionVector = ToEigen(X, dofTypes);
        mImplicitCallBack.SetReducedSolutionSpaceOperator(mReducedSolutionSpaceOperator);
    }

    void Solve(std::string solverString)
    {
        mSolver.DoStep(mSolutionVector, mImplicitCallBack, 0, solverString);
    }

private:
    MeshFem mMesh;
    DofType mDof;

    Laws::LinearElastic<3> mElasticLaw;
    Integrands::MomentumBalance<3> mMomentumBalance;

    IntegrationTypeTensorProduct<3> mIntegration;

    CellStorage mCells;

    TimeDependentProblem mFunctions;
    QuasistaticSolver mSolver;
    ReducedSolutionSpace mReducedSolutionSpaceOperator;
    ImplicitCallBack mImplicitCallBack;

    Eigen::VectorXd mSolutionVector;
};

static void ElasticSolveMumps(benchmark::State& state)
{
    TestStructure s(Eigen::Vector3i(10, 10, 100));
    for (auto _ : state)
        s.Solve("MumpsLDLT");
}

BENCHMARK(ElasticSolveMumps);
BENCHMARK_MAIN();
