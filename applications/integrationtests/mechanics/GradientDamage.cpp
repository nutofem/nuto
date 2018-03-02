#include "BoostUnitTest.h"

#include <boost/filesystem.hpp>

#include "math/EigenCompanion.h"

#include "mechanics/integrands/GradientDamage.h"
#include "mechanics/integrands/NeumannBc.h"
#include "mechanics/constitutive/damageLaws/DamageLawExponential.h"
#include "mechanics/mesh/UnitMeshFem.h"
#include "mechanics/mesh/MeshGmsh.h"
#include "mechanics/mesh/MeshFemDofConvert.h"
#include "mechanics/interpolation/InterpolationTrussLobatto.h"
#include "mechanics/integrationtypes/IntegrationType2D3NGauss12Ip.h"
#include "mechanics/integrationtypes/IntegrationTypeTensorProduct.h"
#include "mechanics/tools/CellStorage.h"
#include "mechanics/tools/TimeDependentProblem.h"
#include "mechanics/tools/QuasistaticSolver.h"
#include "mechanics/tools/AdaptiveSolve.h"
#include "mechanics/constraints/ConstraintCompanion.h"

#include "visualize/PostProcess.h"
#include "visualize/VoronoiGeometries.h"

using namespace NuTo;


BOOST_AUTO_TEST_CASE(Integrand)
{
    DofType d("Displacements", 1);
    ScalarDofType eeq("NonlocalEquivalentStrains");

    /* MATERIAL */
    double E = 30000;
    double nu = 0.2;
    double ft = 4;
    double fc = 40;
    double gf = 0.021;
    double c = 1.25;

    double k0 = ft / E;
    Laws::LinearElastic<1> elasticLaw(E, nu);
    Constitutive::DamageLawExponential dmg(k0, ft / gf, 1.);
    Constitutive::ModifiedMisesStrainNorm<1> strainNorm(nu, fc / ft);

    using Gdm = Integrands::GradientDamage<1, Constitutive::DamageLawExponential>;
    Gdm gdm(d, eeq, c, elasticLaw, dmg, strainNorm);

    /* mesh, interpolations, constraints */
    double L = 40;
    MeshFem mesh = UnitMeshFem::Transform(UnitMeshFem::CreateLines(40),
                                          [&](Eigen::VectorXd x) { return Eigen::VectorXd::Constant(1, x[0] * L); });

    InterpolationTrussLobatto interpolationD(2);
    AddDofInterpolation(&mesh, d, interpolationD);
    AddDofInterpolation(&mesh, eeq);

    Constraint::Constraints constraints;
    constraints.Add(d, Constraint::Component(mesh.NodesAtAxis(eDirection::X, d), {eDirection::X}));
    constraints.Add(d, Constraint::Component(mesh.NodesAtAxis(eDirection::X, d, L), {eDirection::X},
                                             Constraint::RhsRamp(1, 0.2)));

    /* integration cells */
    IntegrationTypeTensorProduct<1> integration(2, eIntegrationMethod::GAUSS);
    const int nIp = integration.GetNumIntegrationPoints();
    CellStorage cellStorage;
    auto cells = cellStorage.AddCells(mesh.ElementsTotal(), integration);

    /* resize history and apply imperfection in the middle */
    gdm.mKappas.resize(cells.Size(), nIp);
    gdm.mKappas.row(cells.Size() / 2) = Eigen::VectorXd::Constant(nIp, 3 * k0);
    gdm.mKappas.row(cells.Size() / 2 + 1) = Eigen::VectorXd::Constant(nIp, 3 * k0);

    /* define time dependent functions */
    TimeDependentProblem equations(&mesh);
    equations.AddGradientFunction(cells, TimeDependentProblem::Bind(gdm, &Gdm::Gradient));
    equations.AddHessian0Function(cells, TimeDependentProblem::Bind(gdm, &Gdm::Hessian0));
    equations.AddUpdateFunction(cells, TimeDependentProblem::Bind(gdm, &Gdm::Update));

    QuasistaticSolver problem(equations, {d, eeq});
    problem.mTolerance = 1.e-6;
    problem.SetConstraints(constraints);

    Visualize::PostProcess visu("GradientDamageOut1D");
    visu.DefineVisualizer("GDM", cells, Visualize::VoronoiHandler(Visualize::VoronoiGeometryLine(2)));
    visu.Add("GDM", d);
    visu.Add("GDM", eeq);
    visu.Add("GDM", [&](const CellIpData& cipd) { return gdm.Kappa(cipd); }, "Kappa");
    visu.Add("GDM", [&](const CellIpData& cipd) { return cipd.Apply(d, Nabla::Strain()); }, "strain");

    /* solve adaptively */
    auto doStep = [&](double t) { return problem.DoStep(t, "MumpsLU"); };
    auto postProcessF = [&](double t) { visu.Plot(t, true); };
    AdaptiveSolve adaptiveSolve(doStep, postProcessF);
    adaptiveSolve.dt = 0.01;
    BOOST_CHECK_NO_THROW(adaptiveSolve.Solve(3.));

    // A small zone around the middle is damaged and the strains localize there.
    // The rest of the structure is expected to be unloaded
    //
    //  damage:
    //           _
    //          ' '
    //         |   |
    //  _______|   |_______
    //
    // displacements:
    //             ________
    //            ,
    //           |
    //           |
    // _________,
    //
    // nonlocal equivalent strains:
    //
    //           |
    //           |
    // _________/ \________
    //
    //
    // ZoneA: from x = [0 .. middle - damaged zone width]
    //      - displacements ~ 0
    //      - eeq ~ 0
    // ZoneB: from x = [middle + damaged zone width .. L]
    //      - displacements ~ 0.6  ( boundary condition 0.2 / s * 3s )
    //      - eeq ~ 0
    //
    auto& dNodeFromZoneA = mesh.NodeAtCoordinate(EigenCompanion::ToEigen(L / 4), d);
    auto& eeqNodeFromZoneA = mesh.NodeAtCoordinate(EigenCompanion::ToEigen(L / 4), eeq);

    BOOST_CHECK_SMALL(dNodeFromZoneA.GetValues()[0], 1.e-4);
    BOOST_CHECK_SMALL(eeqNodeFromZoneA.GetValues()[0], 1.e-4);

    auto& dNodeFromZoneB = mesh.NodeAtCoordinate(EigenCompanion::ToEigen(3 * L / 4), d);
    auto& eeqNodeFromZoneB = mesh.NodeAtCoordinate(EigenCompanion::ToEigen(3 * L / 4), eeq);

    BOOST_CHECK_CLOSE(dNodeFromZoneB.GetValues()[0], 0.6, 1.e-4);
    BOOST_CHECK_SMALL(eeqNodeFromZoneB.GetValues()[0], 1.e-4);
}


BOOST_AUTO_TEST_CASE(Integrand2D)
{
    boost::filesystem::path binaryPath = boost::unit_test::framework::master_test_suite().argv[0];
    binaryPath.remove_filename();
    std::string meshFile = binaryPath.string() + "/meshes/Holes.msh";

    double E = 30000;
    double nu = 0.2;
    Laws::LinearElastic<2> elasticLaw(E, nu);

    double ft = 4;
    double gf = 0.021;
    Constitutive::DamageLawExponential dmg(ft / E, ft / gf, 1.);

    double fc = 40;
    Constitutive::ModifiedMisesStrainNorm<2> strainNorm(nu, fc / ft);


    DofType d("Displacements", 2);
    ScalarDofType eeq("NonlocalEquivalentStrains");

    double c = 0.1;
    using Gdm = Integrands::GradientDamage<2, Constitutive::DamageLawExponential>;
    using Neumann = Integrands::NeumannBc<2>;

    Gdm gdm(d, eeq, c, elasticLaw, dmg, strainNorm);
    Neumann neumann(d, Eigen::Vector2d(.42, .12));

    MeshGmsh gmsh(meshFile);

    auto& mesh = gmsh.GetMeshFEM();
    auto& matrixElements = gmsh.GetPhysicalGroup("matrix");
    auto& leftElements = gmsh.GetPhysicalGroup("left");
    AddDofInterpolation(&mesh, d, matrixElements);
    AddDofInterpolation(&mesh, d, leftElements);
    AddDofInterpolation(&mesh, eeq, matrixElements);
    AddDofInterpolation(&mesh, eeq, leftElements);

    IntegrationType2D3NGauss12Ip integration;
    CellStorage cellStorage;
    auto cells = cellStorage.AddCells(matrixElements, integration);
    auto cellsLeft = cellStorage.AddCells(leftElements, integration);
    gdm.mKappas.setZero(cells.Size(), integration.GetNumIntegrationPoints());

    TimeDependentProblem equations(&mesh);

    equations.AddHessian0Function(cells, TimeDependentProblem::Bind(gdm, &Gdm::Hessian0));
    equations.AddGradientFunction(cells, TimeDependentProblem::Bind(gdm, &Gdm::Gradient));
    equations.AddUpdateFunction(cells, TimeDependentProblem::Bind(gdm, &Gdm::Update));
    equations.AddGradientFunction(cellsLeft, TimeDependentProblem::Bind(neumann, &Neumann::ExternalLoad));

    QuasistaticSolver problem(equations, {d, eeq});

    using namespace Constraint;
    Constraints constraints;
    constraints.Add(d, Component(mesh.NodeAtCoordinate(Eigen::Vector2d::Zero(), d), {eDirection::X}));
    constraints.Add(d, Component(mesh.NodesAtAxis(eDirection::Y, d), {eDirection::Y}));
    constraints.Add(d, Component(mesh.NodesAtAxis(eDirection::Y, d, 16), {eDirection::Y}, RhsRamp(1, 0.01)));

    problem.SetConstraints(constraints);
    problem.mTolerance = 1.e-6;


    using namespace NuTo::Visualize;
    PostProcess visu("./GradientDamageOut2D");
    visu.DefineVisualizer("GDM", cells, VoronoiHandler(VoronoiGeometryTriangle(integration)));
    visu.Add("GDM", d);
    visu.Add("GDM", eeq);
    visu.Add("GDM", [&](const NuTo::CellIpData& data) { return gdm.mDamageLaw.Damage(gdm.Kappa(data)); }, "Damage");

    auto doStep = [&](double t) { return problem.DoStep(t, "MumpsLU"); };
    auto postProcess = [&](double t) { return visu.Plot(t, true); };

    NuTo::AdaptiveSolve adaptive(doStep, postProcess);
    adaptive.dt = 0.01;
    adaptive.Solve(0.01); // Only one step for this test. Increase to 1., if you want to see magic happen.
}
