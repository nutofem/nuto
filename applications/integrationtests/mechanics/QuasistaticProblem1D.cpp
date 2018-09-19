#include "BoostUnitTest.h"

#include "nuto/math/NewtonRaphson.h"
#include "nuto/math/EigenCompanion.h"

#include "nuto/mechanics/dofs/DofNumbering.h"

#include "nuto/mechanics/integrands/MomentumBalance.h"
#include "nuto/mechanics/integrands/Bind.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeTensorProduct.h"
#include "nuto/mechanics/constraints/ConstraintCompanion.h"
#include "nuto/mechanics/mesh/UnitMeshFem.h"
#include "nuto/mechanics/mesh/GeometryMeshFem.h"
#include "nuto/mechanics/mesh/MeshFemDofConvert.h"
#include "nuto/mechanics/tools/CellStorage.h"
#include "nuto/mechanics/tools/QuasistaticSolver.h"
#include "nuto/mechanics/tools/AdaptiveSolve.h"

#include "nuto/mechanics/constitutive/LocalIsotropicDamage.h"
#include "nuto/mechanics/constitutive/damageLaws/DamageLawExponential.h"
#include "nuto/mechanics/constitutive/ModifiedMisesStrainNorm.h"

using namespace NuTo;

class LocalDamageTruss
{
public:
    LocalDamageTruss(int numElements, Material::Softening m)
        : mGeoMesh(UnitMeshFem::CreateLines(numElements))
        , mMesh(mGeoMesh)
        , mDof("Dispacement", 1)
        , mLaw(m)
        , mMomentumBalance(mDof, mLaw)
        , mEquations(&mMesh)
        , mProblem(mEquations, mDof)
        , mIntegrationType(2, eIntegrationMethod::GAUSS)
    {
        AddDofInterpolation(&mMesh, mDof);

        mCellGroup = mCells.AddCells(mMesh.ElementsTotal(), mIntegrationType);
        mLaw.mEvolution.ResizeHistoryData(mCellGroup.Size(), mIntegrationType.GetNumIntegrationPoints());


        auto Gradient = TimeDependentProblem::Bind_dt(mMomentumBalance, &Integrands::MomentumBalance<1>::Gradient);
        auto Hessian0 = TimeDependentProblem::Bind_dt(mMomentumBalance, &Integrands::MomentumBalance<1>::Hessian0);
        TimeDependentProblem::UpdateFunction UpdateHistory = [&](const CellIpData& cellIpData, double, double dt) {
            mLaw.Update(cellIpData.Apply(mDof, Nabla::Strain()), dt, cellIpData.Ids());
        };

        mEquations.AddGradientFunction(mCellGroup, Gradient);
        mEquations.AddHessian0Function(mCellGroup, Hessian0);
        mEquations.AddUpdateFunction(mCellGroup, UpdateHistory);

        auto constraints = DefineConstraints(mMesh, mDof);

        mProblem.SetConstraints(constraints);
        mProblem.mTolerance = 1.e-12;
    }

    void SetImperfection(double kappaImperfection)
    {
        const int imperfectionCell = mMesh.Elements.Size() / 2;
        mLaw.mEvolution.mKappas(imperfectionCell, 0) = kappaImperfection;
    }

    void Solve(double tEnd)
    {
        auto doStep = [&](double t) { return mProblem.DoStep(t, "EigenSparseLU"); };
        AdaptiveSolve adaptive(doStep);
        adaptive.dt = 0.01;
        adaptive.Solve(tEnd);
    }

    std::vector<double> DamageField()
    {
        auto Damage = [&](const CellIpData& cellIpData) {
            double kappa = mLaw.mEvolution.mKappas(cellIpData.Ids().cellId, cellIpData.Ids().ipId);
            return Eigen::VectorXd::Constant(1, mLaw.mDamageLaw.Damage(kappa));
        };

        std::vector<double> damage;
        for (auto& cell : mCellGroup)
        {
            auto ipDamage = cell.Eval(Damage);
            for (auto dmg : ipDamage)
                damage.push_back(dmg[0]);
        }
        return damage;
    }

private:
    GeometryMeshFem mGeoMesh;
    MeshFem mMesh;

    DofType mDof;
    Laws::LocalIsotropicDamage<1> mLaw;
    Integrands::MomentumBalance<1> mMomentumBalance;

    TimeDependentProblem mEquations;
    QuasistaticSolver mProblem;

    IntegrationTypeTensorProduct<1> mIntegrationType;
    CellStorage mCells;
    Group<CellInterface> mCellGroup;

    Constraint::Constraints DefineConstraints(MeshFem& mesh, DofType disp)
    {
        using namespace NuTo::EigenCompanion;
        auto& nodeLeft = mesh.NodeAtCoordinate(ToEigen(0), disp);
        auto& nodeRight = mesh.NodeAtCoordinate(ToEigen(1), disp);

        auto& nodeMiddle = mesh.NodeAtCoordinate(ToEigen(0.4), disp);

        Constraint::Constraints c;
        Constraint::Equation periodic(nodeRight, 0, Constraint::RhsRamp(1., 0.02));
        periodic.AddIndependentTerm({nodeLeft, 0, -1});

        c.Add(disp, Constraint::Component(nodeMiddle, {eDirection::X}));
        c.Add(disp, periodic);
        return c;
    }
};

BOOST_AUTO_TEST_CASE(LocalDamage1DWithoutImperfection)
{
    auto material = Material::DefaultConcrete();
    material.fMin = 0.;
    LocalDamageTruss problem(5, material);
    problem.Solve(1);
    auto damageField = problem.DamageField();
    // without an imperfection, we expect a uniform damage distribution.
    // If the trial state is messed up, we would have localization at
    // the displacement BC. (damage = 1 there, and damage = 0 everywhere else)
    auto firstDamage = damageField[0];
    for (double damage : damageField)
        BOOST_CHECK_CLOSE(damage, firstDamage, 1.e-5);
}

BOOST_AUTO_TEST_CASE(LocalDamage1DWithImperfection)
{
    auto material = Material::DefaultConcrete();
    material.fMin = 0.;
    LocalDamageTruss problem(5, material);
    problem.SetImperfection(0.001);
    problem.Solve(1);
    auto damageField = problem.DamageField();
    // We expect localization in the two integration points
    // near the center
    size_t damagedIp1 = 4;
    size_t damagedIp2 = 5;
    for (size_t i = 0; i < damageField.size(); ++i)
    {
        if (i == damagedIp1 or i == damagedIp2)
            BOOST_CHECK_CLOSE(damageField[i], 1., 1.e-2);
        else
            BOOST_CHECK_SMALL(damageField[i], 1.e-5);
    }
}
