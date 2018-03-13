#include <fstream>
#include <boost/filesystem.hpp>

#include "math/EigenCompanion.h"

#include "mechanics/integrands/GradientDamage.h"
#include "mechanics/integrands/MomentumBalance.h"
#include "mechanics/constitutive/LinearElastic.h"
#include "mechanics/constitutive/LocalIsotropicDamage.h"
#include "mechanics/constitutive/damageLaws/DamageLawExponential.h"
#include "mechanics/mesh/MeshGmsh.h"
#include "mechanics/mesh/MeshFemDofConvert.h"
#include "mechanics/integrationtypes/IntegrationTypeTensorProduct.h"
#include "mechanics/integrationtypes/IntegrationType2D3NGauss6Ip.h"
#include "mechanics/tools/CellStorage.h"
#include "mechanics/tools/TimeDependentProblem.h"
#include "mechanics/tools/QuasistaticSolver.h"
#include "mechanics/tools/AdaptiveSolve.h"
#include "mechanics/constraints/ConstraintCompanion.h"
#include "visualize/PostProcess.h"
#include "visualize/VoronoiGeometries.h"

using namespace NuTo;

int main(int, char* argv[])
{
    boost::filesystem::path binaryPath = argv[0];
    binaryPath.remove_filename();
    std::string meshFile = binaryPath.string() + "/SingleInclusion.msh";

    double E = 26738;
    double nu = 0.18;
    Laws::LinearElasticDamage<2> unilateral(E, nu, false);
    Laws::LinearElastic<2> elasticLaw(2 * E, nu);

    double ft = 3.4;
    double k0 = ft / E;
    double Gf = 0.12; // global
    double gf = 0.0216; // calibrated via 1D tensile
    Constitutive::DamageLawExponential dmg(k0, ft / (gf * 10), 0.99);

    double fc = ft * 10;
    Constitutive::ModifiedMisesStrainNorm<2> strainNorm(nu, fc / ft);

    double F = 0.75;
    double thickness = 0.1;

    double ftItz = ft * F;
    double gfItz = Gf / thickness * F;

    double k0Itz = ftItz / E;
    Constitutive::ModifiedMisesStrainNorm<2> strainNormItz(nu, 1);
    Constitutive::DamageLawExponential iztDmg(k0Itz, ftItz / gfItz, 0.99);
    Laws::LocalIsotropicDamage<2, Constitutive::DamageLawExponential, Laws::EvolutionImplicit<2>> itzLaw(
            unilateral, iztDmg, strainNormItz);

    DofType d("Displacements", 2);
    ScalarDofType eeq("NonlocalEquivalentStrains");

    double c = 2.0;
    using Gdm = Integrands::GradientDamage<2, Constitutive::DamageLawExponential, NonlocalInteraction::Decreasing>;
    Gdm gdm(d, eeq, c, unilateral, dmg, strainNorm, NonlocalInteraction::Decreasing(0.05, 5));
    // using Gdm = Integrands::GradientDamage<2, Constitutive::DamageLawExponential, NonlocalInteraction::Constant>;
    // Gdm gdm(d, eeq, c, unilateral, dmg, strainNorm);

    using Mb = Integrands::MomentumBalance<2>;
    Mb mbAggregates(d, elasticLaw);
    // Mb mbItz(d, elasticLaw);
    Mb mbItz(d, itzLaw);


    MeshGmsh gmsh(meshFile);
    auto& mesh = gmsh.GetMeshFEM();

    AddDofInterpolation(&mesh, d);
    AddDofInterpolation(&mesh, eeq, gmsh.GetPhysicalGroup("Matrix"));

    IntegrationType2D3NGauss6Ip integration;
    IntegrationTypeTensorProduct<2> integrationItz(3, eIntegrationMethod::GAUSS);

    CellStorage cellStorage;
    auto cellsMatrix = cellStorage.AddCells(gmsh.GetPhysicalGroup("Matrix"), integration);
    auto cellsAggregates = cellStorage.AddCells(gmsh.GetPhysicalGroup("Aggregates"), integration);
    auto cellsItz = cellStorage.AddCells(gmsh.GetPhysicalGroup("Interfaces"), integrationItz);

    gdm.mKappas.setZero(cellsMatrix.Size(), integration.GetNumIntegrationPoints());
    itzLaw.mEvolution.mKappas.setZero(cellsItz.Size(), integrationItz.GetNumIntegrationPoints());

    TimeDependentProblem equations(&mesh);

    equations.AddHessian0Function(cellsMatrix, TimeDependentProblem::Bind(gdm, &Gdm::Hessian0));
    equations.AddGradientFunction(cellsMatrix, TimeDependentProblem::Bind(gdm, &Gdm::Gradient));
    equations.AddUpdateFunction(cellsMatrix, TimeDependentProblem::Bind(gdm, &Gdm::Update));

    equations.AddHessian0Function(cellsAggregates, TimeDependentProblem::Bind_dt(mbAggregates, &Mb::Hessian0));
    equations.AddGradientFunction(cellsAggregates, TimeDependentProblem::Bind_dt(mbAggregates, &Mb::Gradient));

    equations.AddHessian0Function(cellsItz, TimeDependentProblem::Bind_dt(mbItz, &Mb::Hessian0));
    equations.AddGradientFunction(cellsItz, TimeDependentProblem::Bind_dt(mbItz, &Mb::Gradient));
    TimeDependentProblem::UpdateFunction UpdateItz = [&](const CellIpData& data, double, double) {
        itzLaw.Update(data.Apply(d, Nabla::Strain()), 0., data.Ids());
    };
    equations.AddUpdateFunction(cellsItz, UpdateItz);

    QuasistaticSolver problem(equations, {d, eeq});

    using namespace Constraint;
    Constraints constraints;
    constraints.Add(d, Component(mesh.NodesAtAxis(eDirection::Y, d), {eDirection::Y}));
    // constraints.Add(d, Component(mesh.NodesAtAxis(eDirection::Y, d), {eDirection::X}));
    constraints.Add(d, Component(mesh.NodeAtCoordinate(Eigen::Vector2d(0, 0), d), {eDirection::X}));
    // constraints.Add(d, Component(mesh.NodeAtCoordinate(Eigen::Vector2d(60, 0), d), {eDirection::X}));

    // constraints.Add(d, Component(mesh.NodeAtCoordinate(Eigen::Vector2d(0, 60), d), {eDirection::X}));
    // constraints.Add(d, Component(mesh.NodeAtCoordinate(Eigen::Vector2d(60, 60), d), {eDirection::X}));

    auto topNodes = mesh.NodesAtAxis(eDirection::Y, d, 60);
    // constraints.Add(d, Component(topNodes, {eDirection::X}, RhsRamp(1, -0.01)));
    constraints.Add(d, Component(topNodes, {eDirection::Y}, RhsRamp(1, -0.1)));

    problem.SetConstraints(constraints);
    problem.mTolerance = 1.e-6;


    using namespace NuTo::Visualize;
    PostProcess visu("./ConcreteCompressionWithInclusionOut");
    visu.DefineVisualizer("Matrix", cellsMatrix, VoronoiHandler(VoronoiGeometryTriangle(integration)));
    visu.Add("Matrix", d);
    visu.Add("Matrix", eeq);
    visu.Add("Matrix", [&](const NuTo::CellIpData& data) { return gdm.mDamageLaw.Damage(gdm.Kappa(data)); }, "Damage");
    visu.Add("Matrix", [&](const NuTo::CellIpData& data) { return gdm.Kappa(data) / k0; }, "k0Factor");


    visu.DefineVisualizer("Itz", cellsItz, VoronoiHandler(VoronoiGeometryQuad(3)));
    visu.Add("Itz", d);
    visu.Add("Itz",
             [&](const NuTo::CellIpData& data) {
                 double kappa = itzLaw.mEvolution.mKappas(data.Ids().cellId, data.Ids().ipId);
                 return itzLaw.mDamageLaw.Damage(kappa);
             },
             "Damage");
    visu.Add("Itz",
             [&](const NuTo::CellIpData& data) {
                 return itzLaw.mEvolution.mKappas(data.Ids().cellId, data.Ids().ipId) / k0Itz;
             },
             "k0Factor");


    std::ofstream loadDisp(visu.ResultDirectory() + "/LD.dat");

    auto doStep = [&](double t) { return problem.DoStep(t, "MumpsLU"); };
    auto postProcess = [&](double t) {
        visu.Plot(t, true);
        problem.WriteTimeDofResidual(loadDisp, d, DofNumbering::Get(topNodes, ToComponentIndex(eDirection::Y)));
    };

    NuTo::AdaptiveSolve adaptive(doStep, postProcess);
    adaptive.dt = 0.1;
    adaptive.Solve(3); // Only one step for this test. Increase to 1., if you want to see magic happen.
}
