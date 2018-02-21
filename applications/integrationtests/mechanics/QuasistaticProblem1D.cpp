#include "BoostUnitTest.h"

#include "math/NewtonRaphson.h"
#include "math/EigenCompanion.h"

#include "mechanics/dofs/DofNumbering.h"

#include "mechanics/integrands/MomentumBalance.h"
#include "mechanics/integrands/Bind.h"
#include "mechanics/integrationtypes/IntegrationTypeTensorProduct.h"
#include "mechanics/constraints/ConstraintCompanion.h"
#include "mechanics/mesh/UnitMeshFem.h"
#include "mechanics/mesh/MeshFemDofConvert.h"
#include "mechanics/tools/CellStorage.h"
#include "mechanics/tools/QuasistaticSolver.h"

#include "mechanics/constitutive/LinearElastic.h"
#include "mechanics/constitutive/LocalIsotropicDamage.h"
#include "mechanics/constitutive/damageLaws/DamageLawExponential.h"
#include "mechanics/constitutive/ModifiedMisesStrainNorm.h"

#include "visualize/Visualizer.h"
#include "visualize/AverageHandler.h"
#include "visualize/AverageGeometries.h"

using namespace NuTo;


class LocalDamageTruss
{
public:
    LocalDamageTruss(int numElements)
        : mMesh(UnitMeshFem::CreateLines(numElements))
        , mDof("Dispacement", 1)
        , mLaw(ConstitutiveLaw())
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
        TimeDependentProblem::UpdateFunction UpdateHistory = [&](const CellData& cellData, const CellIpData& cellIpData,
                                                                 double, double dt) {
            EngineeringStrain<1> strain = cellIpData.B(mDof, B::Strain()) * cellData.GetNodeValues(mDof);
            mLaw.Update(strain, dt, cellData.GetCellId(), cellIpData.GetIpId());
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
        const int imperfectionIp = mLaw.mEvolution.Ip(imperfectionCell, 0);
        mLaw.mEvolution.mKappas[imperfectionIp] = kappaImperfection;
    }

    void Solve(double tEnd)
    {
        double globalTime = 0;
        double timeStep = 0.01;
        double maxTimeStep = 0.1;

        int iStep = 0;
        while (globalTime < tEnd && timeStep > 1.e-6)
        {
            std::cout << "Step " << iStep << " at t = " << globalTime << std::endl;

            auto numIterations = mProblem.DoStep(globalTime + timeStep);

            if (not numIterations)
            {
                timeStep *= 0.5;
                std::cout << "No convergence!\n";
                continue; // without updating the global time
            }

            std::cout << "Converence after " << numIterations.get() << " iterations.\n";
            globalTime += timeStep;
            iStep++;

            if (numIterations < 3)
            {
                timeStep *= 1.5;
                timeStep = std::min(timeStep, maxTimeStep);
                std::cout << "--> Increasing time step to " << timeStep << ".\n";
            }
        }
    }

    std::vector<double> DamageField()
    {
        auto Damage = [&](const CellData& cellData, const CellIpData& cellIpData) {
            int ip = mLaw.mEvolution.Ip(cellData.GetCellId(), cellIpData.GetIpId());
            double kappa = mLaw.mEvolution.mKappas[ip];
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
    using ElasticLaw = Laws::LinearElastic<1>;
    using DamageLaw = Constitutive::DamageLawExponential;
    using StrainNorm = Constitutive::ModifiedMisesStrainNorm<1>;
    using Evolution = Laws::EvolutionImplicit<1, StrainNorm>;
    using Law = Laws::LocalIsotropicDamage<1, DamageLaw, Evolution, ElasticLaw>;

    MeshFem mMesh;

    DofType mDof;
    Law mLaw;
    Integrands::MomentumBalance<1> mMomentumBalance;

    TimeDependentProblem mEquations;
    QuasistaticSolver mProblem;

    IntegrationTypeTensorProduct<1> mIntegrationType;
    CellStorage mCells;
    Group<CellInterface> mCellGroup;

    static Law ConstitutiveLaw()
    {
        double E = 20000.;
        double nu = 0.2;
        double ft = 4;
        double fc = 10 * ft;

        double k0 = ft / E;
        double alpha = 0.99;
        double beta = 350;

        ElasticLaw elasticLaw(E, nu);
        DamageLaw damageLaw(k0, beta, alpha);
        StrainNorm strainNorm(nu, fc / ft);
        Evolution evolution(strainNorm);

        return Law(damageLaw, evolution, elasticLaw);
    }

    Constraint::Constraints DefineConstraints(MeshFem& mesh, DofType disp)
    {
        using namespace NuTo::EigenCompanion;
        auto& nodeLeft = mesh.NodeAtCoordinate(ToEigen(0), disp);
        auto& nodeRight = mesh.NodeAtCoordinate(ToEigen(1), disp);

        auto& nodeMiddle = mesh.NodeAtCoordinate(ToEigen(0.4), disp);

        Constraint::Constraints c;


        Constraint::Equation periodic(nodeRight, 0, Constraint::RhsRamp(1., 0.01));
        periodic.AddTerm({nodeLeft, 0, -1});

        c.Add(disp, Constraint::Component(nodeMiddle, {eDirection::X}));
        c.Add(disp, periodic);
        return c;
    }
};

BOOST_AUTO_TEST_CASE(LocalDamage1DWithoutImperfection)
{
    LocalDamageTruss problem(5);
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
    LocalDamageTruss problem(5);
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
