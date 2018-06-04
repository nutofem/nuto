#include "BoostUnitTest.h"

#include <boost/filesystem.hpp>

#include "nuto/math/EigenCompanion.h"
#include "nuto/mechanics/dofs/DofVectorConvertEigen.h"

#include "nuto/mechanics/integrands/GradientDamage.h"
#include "nuto/mechanics/integrands/NeumannBc.h"
#include "nuto/mechanics/mesh/UnitMeshFem.h"
#include "nuto/mechanics/mesh/MeshGmsh.h"
#include "nuto/mechanics/mesh/MeshFemDofConvert.h"
#include "nuto/mechanics/interpolation/InterpolationTrussLobatto.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeTriangle.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeTensorProduct.h"
#include "nuto/mechanics/tools/CellStorage.h"
#include "nuto/mechanics/tools/TimeDependentProblem.h"
#include "nuto/mechanics/tools/QuasistaticSolver.h"
#include "nuto/mechanics/tools/AdaptiveSolve.h"
#include "nuto/mechanics/constraints/ConstraintCompanion.h"
#include "nuto/visualize/PostProcess.h"
#include "nuto/visualize/VoronoiGeometries.h"

using namespace NuTo;

BOOST_AUTO_TEST_CASE(Integrand)
{
    DofType d("Displacements", 1);
    ScalarDofType eeq("NonlocalEquivalentStrains");

    double L = 40;
    auto material = Material::DefaultConcrete();
    using Gdm = Integrands::GradientDamage<1, NonlocalInteraction::Decreasing>;
    // The global fracture energy Gf is influenced by the nonlocal parameter. A smaller nonlocal parameter (due to
    // decreasing interaction) results in a smaller Gf. Not adapting the material parameter gf will result in
    //  1) an unexpected (wrong) global Gf (not crucial for this test)
    //  2) a snap-back that causes the direct displacement controlled test to not converge. Problem.
    // The factor 10 is by no means accurate. It just avoids the snap back.
    double gfNonlocalDecreasingInteractionFactor = 10;
    material.gf *= gfNonlocalDecreasingInteractionFactor;
    material.c = 0.25;
    double k0 = material.ft / material.E;

    Gdm gdm(d, eeq, material);

    /* mesh, interpolations, constraints */
    MeshFem mesh = UnitMeshFem::Transform(UnitMeshFem::CreateLines(80),
                                          [&](Eigen::VectorXd x) { return Eigen::VectorXd::Constant(1, x[0] * L); });

    InterpolationTrussLobatto interpolationD(2);
    AddDofInterpolation(&mesh, d, interpolationD);
    AddDofInterpolation(&mesh, eeq);

    Constraint::Constraints constraints;
    constraints.Add(d, Constraint::Component(mesh.NodesAtAxis(eDirection::X, d), {eDirection::X}));
    constraints.Add(d, Constraint::Component(mesh.NodesAtAxis(eDirection::X, d, L), {eDirection::X},
                                             Constraint::RhsRamp(1, 0.2)));

    /* integration cells */
    IntegrationTypeTensorProduct<1> integration(3, eIntegrationMethod::GAUSS);
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

    std::vector<DofType> dofTypes;
    dofTypes.push_back(d);
    dofTypes.push_back(eeq);

    DofVector<double> X = equations.RenumberDofs(constraints, dofTypes, DofVector<double>());

    DofContainer<int> numTotalDofs;
    DofInfo dofInfoDisp = DofNumbering::Build(mesh.NodesTotal(d), d, constraints);
    DofInfo dofInfoEeq = DofNumbering::Build(mesh.NodesTotal(eeq), eeq, constraints);
    numTotalDofs.Insert(d, dofInfoDisp.numDependentDofs[d] + dofInfoDisp.numIndependentDofs[d]);
    numTotalDofs.Insert(eeq, dofInfoEeq.numDependentDofs[eeq] + dofInfoEeq.numIndependentDofs[eeq]);

    ReducedSolutionSpace reducedSolutionSpaceOperator(dofTypes, numTotalDofs, constraints);

    ImplicitCallBack implicitCallBack(equations, reducedSolutionSpaceOperator);

    QuasistaticSolver problem(X);

    int dofLeft = mesh.NodeAtCoordinate(EigenCompanion::ToEigen(L), d).GetDofNumber(0);

    Visualize::PostProcess visu("GradientDamageOut1D");
    visu.DefineVisualizer("GDM", cells, Visualize::VoronoiHandler(Visualize::VoronoiGeometryLine(3)));
    visu.Add("GDM", d);
    visu.Add("GDM", eeq);
    visu.Add("GDM", [&](const CellIpData& cipd) { return gdm.Kappa(cipd); }, "Kappa");
    visu.Add("GDM", [&](const CellIpData& cipd) { return gdm.mDamageLaw.Damage(gdm.Kappa(cipd)); }, "Damage");
    visu.Add("GDM", [&](const CellIpData& cipd) { return cipd.Apply(d, Nabla::Strain()); }, "strain");
    visu.Add("GDM",
             [&](const CellIpData& cipd) {
                 double omega = gdm.mDamageLaw.Damage(gdm.Kappa(cipd));
                 double R = 0.005;
                 double N = 5;
                 return ((1. - R) * std::exp(-N * omega) + R - std::exp(-N)) / (1 - std::exp(-N));
             },
             "g");

    std::ofstream loadDisplacement(visu.ResultDirectory() + "/LD.dat");

    /* solve adaptively */
    auto doStep = [&](double t) { return problem.DoStep(implicitCallBack, t, "MumpsLU", 1.e-6); };
    auto postProcessF = [&](double t) {
        visu.Plot(t, true);
        problem.WriteTimeDofResidual(loadDisplacement, d, {dofLeft}, implicitCallBack);
    };

    AdaptiveSolve adaptiveSolve(doStep, postProcessF);
    adaptiveSolve.dt = 0.01;
    adaptiveSolve.Solve(3.);

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
    auto& dNodeFromZoneA = mesh.NodeAtCoordinate(EigenCompanion::ToEigen(L / 8), d);
    auto& eeqNodeFromZoneA = mesh.NodeAtCoordinate(EigenCompanion::ToEigen(L / 8), eeq);

    BOOST_CHECK_SMALL(dNodeFromZoneA.GetValues()[0], 1.e-4);
    BOOST_CHECK_SMALL(eeqNodeFromZoneA.GetValues()[0], 1.e-4);

    auto& dNodeFromZoneB = mesh.NodeAtCoordinate(EigenCompanion::ToEigen(7 * L / 8), d);
    auto& eeqNodeFromZoneB = mesh.NodeAtCoordinate(EigenCompanion::ToEigen(7 * L / 8), eeq);

    BOOST_CHECK_CLOSE_FRACTION(dNodeFromZoneB.GetValues()[0], 0.6, 1.e-3);
    BOOST_CHECK_SMALL(eeqNodeFromZoneB.GetValues()[0], 1.e-4);
}


BOOST_AUTO_TEST_CASE(Integrand2D)
{
    boost::filesystem::path binaryPath = boost::unit_test::framework::master_test_suite().argv[0];
    binaryPath.remove_filename();
    std::string meshFile = binaryPath.string() + "/meshes/Holes.msh";

    auto material = Material::DefaultConcrete();
    material.c = 2;

    DofType d("Displacements", 2);
    ScalarDofType eeq("NonlocalEquivalentStrains");

    using Gdm = Integrands::GradientDamage<2>;
    using Neumann = Integrands::NeumannBc<2>;

    Gdm gdm(d, eeq, material);
    Neumann neumann(d, Eigen::Vector2d(.42, .12));

    MeshGmsh gmsh(meshFile);

    auto& mesh = gmsh.GetMeshFEM();
    auto& matrixElements = gmsh.GetPhysicalGroup("matrix");
    auto& leftElements = gmsh.GetPhysicalGroup("left");
    AddDofInterpolation(&mesh, d, matrixElements);
    AddDofInterpolation(&mesh, d, leftElements);
    AddDofInterpolation(&mesh, eeq, matrixElements);
    AddDofInterpolation(&mesh, eeq, leftElements);

    IntegrationTypeTriangle triangleIntegration(4);
    IntegrationTypeTensorProduct<1> lineIntegration(2, eIntegrationMethod::GAUSS);
    CellStorage cellStorage;
    auto cells = cellStorage.AddCells(matrixElements, triangleIntegration);
    auto cellsLeft = cellStorage.AddCells(leftElements, lineIntegration);
    gdm.mKappas.setZero(cells.Size(), triangleIntegration.GetNumIntegrationPoints());

    TimeDependentProblem equations(&mesh);

    equations.AddHessian0Function(cells, TimeDependentProblem::Bind(gdm, &Gdm::Hessian0));
    equations.AddGradientFunction(cells, TimeDependentProblem::Bind(gdm, &Gdm::Gradient));
    equations.AddUpdateFunction(cells, TimeDependentProblem::Bind(gdm, &Gdm::Update));
    equations.AddGradientFunction(cellsLeft, TimeDependentProblem::Bind(neumann, &Neumann::ExternalLoad));

    using namespace Constraint;
    Constraints constraints;
    constraints.Add(d, Component(mesh.NodeAtCoordinate(Eigen::Vector2d::Zero(), d), {eDirection::X}));
    constraints.Add(d, Component(mesh.NodesAtAxis(eDirection::Y, d), {eDirection::Y}));
    auto topNodes = mesh.NodesAtAxis(eDirection::Y, d, 16);
    constraints.Add(d, Component(topNodes, {eDirection::Y}, RhsRamp(1, 0.01)));

    std::vector<DofType> dofTypes;
    dofTypes.push_back(d);
    dofTypes.push_back(eeq);

    DofVector<double> X = equations.RenumberDofs(constraints, dofTypes, DofVector<double>());
    QuasistaticSolver problem(X);

    DofContainer<int> numTotalDofs;
    DofInfo dofInfoDisp = DofNumbering::Build(mesh.NodesTotal(d), d, constraints);
    DofInfo dofInfoEeq = DofNumbering::Build(mesh.NodesTotal(eeq), eeq, constraints);
    numTotalDofs.Insert(d, dofInfoDisp.numDependentDofs[d] + dofInfoDisp.numIndependentDofs[d]);
    numTotalDofs.Insert(eeq, dofInfoEeq.numDependentDofs[eeq] + dofInfoEeq.numIndependentDofs[eeq]);

    ReducedSolutionSpace reducedSolutionSpaceOperator(dofTypes, numTotalDofs, constraints);

    ImplicitCallBack implicitCallBack(equations, reducedSolutionSpaceOperator);

    using namespace NuTo::Visualize;
    PostProcess visu("./GradientDamageOut2D");
    visu.DefineVisualizer("GDM", cells, VoronoiHandler(VoronoiGeometryTriangle(triangleIntegration)));
    visu.Add("GDM", d);
    visu.Add("GDM", eeq);
    visu.Add("GDM", [&](const NuTo::CellIpData& data) { return gdm.mDamageLaw.Damage(gdm.Kappa(data)); }, "Damage");

    std::ofstream loadDisp(visu.ResultDirectory() + "/LD.dat");

    auto doStep = [&](double t) { return problem.DoStep(implicitCallBack, t, "EigenSparseLU", 1.e-6); };
    auto postProcess = [&](double t) {
        visu.Plot(t, true);
        problem.WriteTimeDofResidual(loadDisp, d, DofNumbering::Get(topNodes, ToComponentIndex(eDirection::Y)),
                                     implicitCallBack);
    };

    NuTo::AdaptiveSolve adaptive(doStep, postProcess);
    adaptive.dt = 0.01;
    adaptive.Solve(adaptive.dt); // Only one step for this test. Increase to 1., if you want to see magic happen.
}
