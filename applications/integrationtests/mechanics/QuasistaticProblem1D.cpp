#include "BoostUnitTest.h"

#include "math/NewtonRaphson.h"

#include "mechanics/dofs/DofNumbering.h"

#include "mechanics/integrands/MomentumBalance.h"
#include "mechanics/integrands/Bind.h"
#include "mechanics/integrationtypes/IntegrationTypeTensorProduct.h"
#include "mechanics/constraints/ConstraintCompanion.h"
#include "mechanics/mesh/UnitMeshFem.h"
#include "mechanics/mesh/MeshFemDofConvert.h"
#include "mechanics/tools/CellStorage.h"
#include "mechanics/tools/QuasistaticProblem.h"

#include "mechanics/constitutive/LinearElastic.h"
#include "mechanics/constitutive/LocalIsotropicDamage.h"
#include "mechanics/constitutive/damageLaws/DamageLawExponential.h"
#include "mechanics/constitutive/ModifiedMisesStrainNorm.h"

#include "visualize/Visualizer.h"
#include "visualize/AverageHandler.h"
#include "visualize/AverageGeometries.h"

using namespace NuTo;

auto ConstitutiveLaw(int numCells, int numIpsPerCell)
{
    using ElasticLaw = Laws::LinearElastic<1>;
    using DamageLaw = Constitutive::DamageLawExponential;
    using StrainNorm = Constitutive::ModifiedMisesStrainNorm<1>;
    using Evolution = Laws::EvolutionImplicit<1, StrainNorm>;
    using Law = Laws::LocalIsotropicDamage<1, DamageLaw, Evolution, ElasticLaw>;

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
    Evolution evolution(strainNorm, numCells, numIpsPerCell);

    return Law(damageLaw, evolution, elasticLaw);
}

Constraint::Constraints DefineConstraints(MeshFem& mesh, DofType disp)
{
    auto nodesLeft = mesh.NodesAtAxis(eDirection::X, disp);
    auto nodesRight = mesh.NodesAtAxis(eDirection::X, disp, 1);

    Constraint::Constraints c;
    c.Add(disp, Constraint::Component(nodesLeft, {eDirection::X}));
    c.Add(disp, Constraint::Component(nodesRight, {eDirection::X}, Constraint::RhsRamp(1., 0.01)));
    return c;
}

std::vector<double> DamageField1D(int numElements, double initialKappa)
{
    MeshFem mesh = UnitMeshFem::CreateLines(numElements);
    Group<ElementCollectionFem> elements;
    for (auto& element : mesh.Elements)
        elements.Add(element);

    DofType disp("Displacement", 1);
    AddDofInterpolation(&mesh, disp);

    auto constraints = DefineConstraints(mesh, disp);
    auto dofInfo = DofNumbering::Build(mesh.NodesTotal(disp), disp, constraints);

    CellStorage cells;
    IntegrationTypeTensorProduct<1> integrationType(2, eIntegrationMethod::GAUSS);

    auto cellGroup = cells.AddCells(elements, integrationType);
    auto law = ConstitutiveLaw(cellGroup.Size(), integrationType.GetNumIntegrationPoints());

    // impose imperfection:
    const int imperfectionCell = numElements / 2;
    const int imperfectionIp = law.mEvolution.Ip(imperfectionCell, 0);
    law.mEvolution.mKappas[imperfectionIp] = initialKappa;

    Integrands::MomentumBalance<1> momenumBalance(disp, law);

    auto Gradient = Bind(momenumBalance, &Integrands::MomentumBalance<1>::Gradient);
    auto Hessian0 = Bind(momenumBalance, &Integrands::MomentumBalance<1>::Hessian0);

    CellInterface::VoidFunction UpdateHistory = [&](const CellData& cellData, const CellIpData& cellIpData) {
        EngineeringStrain<1> strain = cellIpData.GetBMatrixStrain(disp) * cellData.GetNodeValues(disp);
        law.Update(strain, 0, cellData.GetCellId(), cellIpData.GetIpId());
    };

    SimpleAssembler assembler(dofInfo.numIndependentDofs, dofInfo.numDependentDofs);
    Merger merger(&mesh);
    EquationSystem system(&assembler, &merger);
    system.AddGradientFunction(cellGroup, Gradient);
    system.AddHessian0Function(cellGroup, Hessian0);
    system.AddUpdateFunction(cellGroup, UpdateHistory);

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> pureSolver;
    auto solver = NewtonRaphson::CreateWrappedSolver(pureSolver);
    Eigen::VectorXd x = Eigen::VectorXd::Zero(dofInfo.numIndependentDofs[disp]);

    double globalTime = 0;
    double timeStep = 0.01;
    double maxTimeStep = 0.1;

    QuasistaticProblem problem(system, constraints, dofInfo.numIndependentDofs[disp], disp);
    problem.mTolerance = 1.e-12;

    int iStep = 0;
    while (globalTime < 1 && timeStep > 1.e-6)
    {
        std::cout << "Step " << iStep << " at t = " << globalTime << std::endl;
        try
        {
            auto trialSystem = problem.TrialSystem(x, globalTime, timeStep);
            Eigen::VectorXd trialX = x + solver.Solve(trialSystem.first, trialSystem.second);

            problem.SetGlobalTime(globalTime + timeStep);

            int numIterations = 0;
            // Some solvers return values > 1e30 as valid solutions. Thus, we store the
            // solution into tmpX (instead of x directly) and check afterwards.
            Eigen::VectorXd tmpX =
                    NewtonRaphson::Solve(problem, trialX, solver, 20., NewtonRaphson::LineSearch(), &numIterations);
            if (tmpX.norm() > 1.e30)
                throw NewtonRaphson::NoConvergence("", "floating point exception");

            problem.UpdateHistory(tmpX);

            globalTime += timeStep;
            x = tmpX;
            iStep++;

            std::cout << "Convergence after " << numIterations << " iterations.\n";
            if (numIterations < 3)
            {
                timeStep *= 1.5;
                timeStep = std::min(timeStep, maxTimeStep);
                std::cout << "--> Increasing time step to " << timeStep << ".\n";
            }
        }
        catch (NewtonRaphson::NoConvergence& e)
        {
            timeStep *= 0.5;
            std::cout << "No convergence!\n";
            std::cout << "--> Reducing time step to " << timeStep << ".\n";
        }
    }

    auto Damage = [&](const CellData& cellData, const CellIpData& cellIpData) {
        int ip = law.mEvolution.Ip(cellData.GetCellId(), cellIpData.GetIpId());
        double kappa = law.mEvolution.mKappas[ip];
        return Eigen::VectorXd::Constant(1, law.mDamageLaw.Damage(kappa));
    };
    std::vector<double> damage;
    for (auto& cell : cellGroup)
    {
        auto ipDamage = cell.Eval(Damage);
        for (auto dmg : ipDamage)
            damage.push_back(dmg[0]);
    }

    return damage;
}

BOOST_AUTO_TEST_CASE(LocalDamage1DWithoutImperfection)
{
    auto damageField = DamageField1D(5., 0.);
    // without an imperfection, we expect a uniform damage distribution.
    // If the trial state is messed up, we would have localization at
    // the displacement BC. (damage = 1 there, and damage = 0 everywhere else)
    auto firstDamage = damageField[0];
    for (double damage : damageField)
        BOOST_CHECK_CLOSE(damage, firstDamage, 1.e-5);
}

BOOST_AUTO_TEST_CASE(LocalDamage1DWithImperfection)
{
    auto damageField = DamageField1D(5., 0.001);
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
