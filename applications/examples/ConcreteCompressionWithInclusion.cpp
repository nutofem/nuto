#include <fstream>
#include <boost/filesystem.hpp>

#include "math/EigenCompanion.h"

#include "mechanics/integrands/GradientDamage.h"
#include "mechanics/integrands/MomentumBalance.h"
#include "mechanics/constitutive/LinearElastic.h"
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

    double E = 30000;
    double nu = 0.2;
    Laws::LinearElasticDamage<2> unilateral(E, nu, true);
    Laws::LinearElastic<2> elasticLaw(2 * E, nu);

    double ft = 4;
    double gf = 0.021 * 10;
    Constitutive::DamageLawExponential dmg(ft / E, ft / gf, 0.99);

    double fc = 4;
    Constitutive::ModifiedMisesStrainNorm<2> strainNorm(nu, fc / ft);


    DofType d("Displacements", 2);
    ScalarDofType eeq("NonlocalEquivalentStrains");

    double c = 1.0;
    using Gdm = Integrands::GradientDamage<2, Constitutive::DamageLawExponential, Integrands::DecreasingInteraction>;
    Gdm gdm(d, eeq, c, unilateral, dmg, strainNorm);

    using Mb = Integrands::MomentumBalance<2>;
    Mb mb(d, elasticLaw);


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

    TimeDependentProblem equations(&mesh);

    equations.AddHessian0Function(cellsMatrix, TimeDependentProblem::Bind(gdm, &Gdm::Hessian0));
    equations.AddGradientFunction(cellsMatrix, TimeDependentProblem::Bind(gdm, &Gdm::Gradient));
    equations.AddUpdateFunction(cellsMatrix, TimeDependentProblem::Bind(gdm, &Gdm::Update));

    equations.AddHessian0Function(cellsAggregates, TimeDependentProblem::Bind_dt(mb, &Mb::Hessian0));
    equations.AddGradientFunction(cellsAggregates, TimeDependentProblem::Bind_dt(mb, &Mb::Gradient));

    equations.AddHessian0Function(cellsItz, TimeDependentProblem::Bind_dt(mb, &Mb::Hessian0));
    equations.AddGradientFunction(cellsItz, TimeDependentProblem::Bind_dt(mb, &Mb::Gradient));

    QuasistaticSolver problem(equations, {d, eeq});

    using namespace Constraint;
    Constraints constraints;
    constraints.Add(d, Component(mesh.NodesAtAxis(eDirection::Y, d), {eDirection::X, eDirection::Y}));

    auto topNodes = mesh.NodesAtAxis(eDirection::Y, d, 60);
    // constraints.Add(d, Component(topNodes, {eDirection::X}, RhsRamp(1, -0.01)));
    constraints.Add(d, Component(topNodes, {eDirection::Y}, RhsRamp(1, -0.1)));

    problem.SetConstraints(constraints);
    problem.mTolerance = 1.e-6;


    using namespace NuTo::Visualize;
    PostProcess visu("./ConcreteCompressionWithInclusionOut");
    visu.DefineVisualizer("GDM", cellsMatrix, VoronoiHandler(VoronoiGeometryTriangle(integration)));
    visu.Add("GDM", d);
    visu.Add("GDM", eeq);
    visu.Add("GDM", [&](const NuTo::CellIpData& data) { return gdm.mDamageLaw.Damage(gdm.Kappa(data)); }, "Damage");

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
