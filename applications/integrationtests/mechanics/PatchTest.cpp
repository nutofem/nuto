#include "BoostUnitTest.h"

#include "boost/ptr_container/ptr_vector.hpp"

#include "nuto/math/EigenCompanion.h"
#include "nuto/math/EigenSparseSolve.h"

#include "nuto/base/Group.h"

#include "nuto/mechanics/dofs/DofNumbering.h"

#include "nuto/mechanics/mesh/MeshFem.h"
#include "nuto/mechanics/mesh/GeometryMeshFem.h"
#include "nuto/mechanics/mesh/MeshFemDofConvert.h"

#include "nuto/mechanics/integrationtypes/IntegrationTypeTensorProduct.h"
#include "nuto/mechanics/interpolation/InterpolationQuadLinear.h"
#include "nuto/mechanics/interpolation/InterpolationTrussLinear.h"

#include "nuto/mechanics/constraints/ConstraintCompanion.h"
#include "nuto/mechanics/constraints/Constraints.h"

#include "nuto/mechanics/constitutive/LinearElastic.h"

#include "nuto/mechanics/integrands/Bind.h"
#include "nuto/mechanics/integrands/MomentumBalance.h"
#include "nuto/mechanics/integrands/NeumannBc.h"

#include "nuto/mechanics/cell/Cell.h"
#include "nuto/mechanics/cell/SimpleAssembler.h"

#include "nuto/visualize/AverageHandler.h"
#include "nuto/visualize/PointHandler.h"
#include "nuto/visualize/Visualizer.h"
#include "nuto/visualize/VoronoiGeometries.h"

#include "nuto/mechanics/solver/Solve.h"

using namespace NuTo;

GeometryMeshFem QuadPatchTestMesh()
{
    /* Something like this:
     *
     *    3-----------------------2
     * /| | - _        e2       / | -->
     * /| |     -7------------6   | -->
     * /| | e4  /     e3      |   | -->
     * /| |    /              |e1 | -->
     * /| |   /    _____------5   | -->
     * /| |  4-----            \  | --> p
     * /| | /      e0           \ | -->
     * /| |/                     \| -->
     *    0-----------------------1
     *   /_\
     *   ///
     *              (c) ttitsche :)
     */
    GeometryMeshFem geoMesh;

    auto& n0 = geoMesh.AddNode(Eigen::Vector2d(0, 0));
    auto& n1 = geoMesh.AddNode(Eigen::Vector2d(10, 0));
    auto& n2 = geoMesh.AddNode(Eigen::Vector2d(10, 10));
    auto& n3 = geoMesh.AddNode(Eigen::Vector2d(0, 10));

    auto& n4 = geoMesh.AddNode(Eigen::Vector2d(2, 2));
    auto& n5 = geoMesh.AddNode(Eigen::Vector2d(8, 3));
    auto& n6 = geoMesh.AddNode(Eigen::Vector2d(8, 7));
    auto& n7 = geoMesh.AddNode(Eigen::Vector2d(4, 7));

    const InterpolationSimple& interpolation = geoMesh.CreateInterpolation(InterpolationQuadLinear());

    geoMesh.AddElement({n0, n1, n5, n4}, interpolation);
    geoMesh.AddElement({n1, n2, n6, n5}, interpolation);
    geoMesh.AddElement({n7, n6, n2, n3}, interpolation);
    geoMesh.AddElement({n4, n5, n6, n7}, interpolation);
    geoMesh.AddElement({n0, n4, n7, n3}, interpolation);

    geoMesh.CreateInterpolation(InterpolationQuadLinear());

    return geoMesh;
}

Constraint::Constraints DefineConstraints(MeshFem* rMesh, DofType dof)
{
    Constraint::Constraints constraints;

    Group<DofNode> nodesConstrainedInX = rMesh->NodesAtAxis(eDirection::X, dof);
    Group<DofNode> nodesConstrainedInY = Group<DofNode>(rMesh->NodeAtCoordinate(Eigen::Vector2d(0, 0), dof));

    constraints.Add(dof, Constraint::Component(nodesConstrainedInX, {eDirection::X}));
    constraints.Add(dof, Constraint::Component(nodesConstrainedInY, {eDirection::Y}));

    return constraints;
}

BOOST_AUTO_TEST_CASE(PatchTestForce)
{
    GeometryMeshFem geoMesh = QuadPatchTestMesh();

    // manually add the boundary element
    const InterpolationSimple& interpolationBc = geoMesh.CreateInterpolation(InterpolationTrussLinear());
    geoMesh.CreateInterpolation(InterpolationQuadLinear());

    // extract existing nodes
    Group<const CoordinateNode> boundaryCoordNodes = geoMesh.NodesAtAxis(eDirection::X, 10);
    const CoordinateNode& nc1 = *boundaryCoordNodes.begin();
    const CoordinateNode& nc2 = *(boundaryCoordNodes.begin() + 1);

    // add the boundary element
    auto& cElmB = geoMesh.AddElement({nc1, nc2}, interpolationBc);

    MeshFem mesh(geoMesh);

    DofType displ("displacements", 2);

    AddDofInterpolation(&mesh, displ);

    Constraint::Constraints constraints = DefineConstraints(&mesh, displ);
    DofInfo dofInfo = DofNumbering::Build(mesh.NodesTotal(displ), displ, constraints);

    // ************************************************************************
    //                 add continuum cells
    // ************************************************************************
    constexpr double E = 20000;
    constexpr double nu = 0.2;
    Laws::LinearElastic<2> linearElasticLaw(E, nu, ePlaneState::PLANE_STRESS);
    Integrands::MomentumBalance<2> momentumBalance(displ, linearElasticLaw);
    auto MomentumGradientF = Bind(momentumBalance, &Integrands::MomentumBalance<2>::Gradient);
    auto MomentumHessian0F = Bind(momentumBalance, &Integrands::MomentumBalance<2>::Hessian0);

    boost::ptr_vector<CellInterface> cellContainer;
    IntegrationTypeTensorProduct<2> integrationType(2, eIntegrationMethod::GAUSS);

    Group<CellInterface> momentumBalanceCells;
    int cellId = 0;
    for (size_t i = 0; i < mesh.NumElements(); i++)
    {
        ElementCollection& element = mesh.GetElement(i);
        if (element.GetShape().Enum() == NuTo::eShape::Quadrilateral)
        {
            cellContainer.push_back(new Cell(element, integrationType, cellId++));
            momentumBalanceCells.Add(cellContainer.back());
        }
    }

    // ************************************************************************
    //                 add boundary cells
    // ************************************************************************

    IntegrationTypeTensorProduct<1> integrationTypeBc(1, eIntegrationMethod::GAUSS);
    Eigen::Vector2d pressureBC(1, 0);
    Integrands::NeumannBc<2> neumannBc(displ, pressureBC);
    auto NeumannLoad = Bind(neumannBc, &Integrands::NeumannBc<2>::ExternalLoad);

    cellContainer.push_back(new Cell(mesh.GetElement(cElmB.GetId()), integrationTypeBc, cellId++));
    auto& neumannCell = cellContainer.back();

    // ************************************************************************
    //                  assemble and solve
    // ************************************************************************
    SimpleAssembler assembler(dofInfo);

    DofVector<double> gradient = assembler.BuildVector(momentumBalanceCells, {displ}, MomentumGradientF);
    // the sign should change here and in the integrand Neumann
    gradient -= assembler.BuildVector({neumannCell}, {displ}, NeumannLoad);

    DofMatrixSparse<double> hessian = assembler.BuildMatrix(momentumBalanceCells, {displ}, MomentumHessian0F);
    DofVector<double> solution = Solve(hessian, -1.0 * gradient, constraints, {displ});

    // Merge dof values %%%%%%%%%%%%%%%%%

    for (DofNode& node : mesh.NodesTotal(displ))
    {
        int dofNumber = node.GetDofNumber(0);
        node.SetValue(0, solution[displ][dofNumber]);
        dofNumber = node.GetDofNumber(1);
        node.SetValue(1, solution[displ][dofNumber]);
    }

    int pointsPerDirection = std::lround(std::sqrt(integrationTypeBc.GetNumIntegrationPoints()));
    pointsPerDirection += 1; // one point per direction doesn't do much Voronoiying
    Visualize::Visualizer visualize(momentumBalanceCells,
                                    Visualize::VoronoiHandler(Visualize::VoronoiGeometryQuad(pointsPerDirection)));
    visualize.DofValues(displ);

    auto stress = [linearElasticLaw, displ](const CellIpData& cellIpData) {
        return linearElasticLaw.Stress(cellIpData.Apply(displ, Nabla::Strain()));
    };
    visualize.CellData(stress, "Stress");

    visualize.CellData([](const CellIpData&) { return EigenCompanion::ToEigen(7.0); }, "Seven");
    visualize.CellData([](const CellIpData& cipd) { return EigenCompanion::ToEigen(cipd.Ids().cellId); }, "CellId");
    visualize.WriteVtuFile("outputVoronoi.vtu");

    auto analyticDisplacementField = [=](Eigen::Vector2d coord) {
        // pressure / A = sigma = E * strain = E * deltaU / L
        return Eigen::Vector2d(pressureBC[0] / E * coord[0], -nu * pressureBC[0] / E * coord[1]);
    };

    for (const CoordinateNode& node : geoMesh.NodesTotal())
    {
        Eigen::VectorXd coord = node.GetCoordinates();
        DofNode& displNode = mesh.NodeAtCoordinate(coord, displ);

        Eigen::VectorXd analyticSolution = analyticDisplacementField(coord);

        BOOST_TEST_MESSAGE("Node at " << coord.transpose() << " with dofs " << displNode.GetDofNumber(0) << ","
                                      << displNode.GetDofNumber(1));
        BOOST_TEST_MESSAGE("( " << displNode.GetValues().transpose() << " )");

        BoostUnitTest::CheckEigenMatrix(displNode.GetValues(), analyticSolution);
    }
}

BOOST_AUTO_TEST_CASE(PatchTestDispl)
{
    GeometryMeshFem geoMesh = QuadPatchTestMesh();
    MeshFem mesh(geoMesh);
    DofType displ("displacements", 2);
    const auto& interpolation = mesh.CreateInterpolation(InterpolationQuadLinear());

    AddDofInterpolation(&mesh, displ, interpolation);

    auto constraints = DefineConstraints(&mesh, displ); // fixed boundary conditions
    Group<DofNode> rightBoundary = mesh.NodesAtAxis(eDirection::X, displ, 10);
    const double boundaryDisplacement = 1.;
    constraints.Add(displ, Constraint::Component(rightBoundary, {eDirection::X}, boundaryDisplacement));

    DofInfo dofInfo = DofNumbering::Build(mesh.NodesTotal(displ), displ, constraints);
    const int numDofs = dofInfo.numIndependentDofs[displ] + dofInfo.numDependentDofs[displ];
    Eigen::SparseMatrix<double> cMatUnit = constraints.BuildUnitConstraintMatrix(displ, numDofs);

    Eigen::MatrixXd cMatUnitFull(cMatUnit);
    BOOST_TEST_MESSAGE("cMatUnit \n" << cMatUnitFull);

    // ************************************************************************
    //   add continuum cells - TODO function to create cells
    // ************************************************************************
    constexpr double E = 20000;
    constexpr double nu = 0.0;
    Laws::LinearElastic<2> linearElasticLaw(E, nu, ePlaneState::PLANE_STRESS);
    Integrands::MomentumBalance<2> momentumBalance(displ, linearElasticLaw);
    auto GradientF = Bind(momentumBalance, &Integrands::MomentumBalance<2>::Gradient);
    auto Hessian0F = Bind(momentumBalance, &Integrands::MomentumBalance<2>::Hessian0);

    boost::ptr_vector<CellInterface> cellContainer;
    IntegrationTypeTensorProduct<2> integrationType(2, eIntegrationMethod::GAUSS);

    Group<CellInterface> cellGroup;
    int cellId = 0;
    for (size_t i = 0; i < mesh.NumElements(); i++)
    {
        auto& element = mesh.GetElement(i);
        cellContainer.push_back(new Cell(element, integrationType, cellId++));
        cellGroup.Add(cellContainer.back());
    }

    // *************************************************************************************
    //   compute initial state (for nonlinear problems, do a trial state
    //   computation first)
    // *************************************************************************************
    DofVector<double> u;
    double t(0);
    u[displ] = constraints.GetSparseGlobalRhs(displ, numDofs, t);
    for (auto& node : mesh.NodesTotal(displ))
    {
        node.SetValue(0, u(displ, node.GetDofNumber(0)));
        node.SetValue(1, u(displ, node.GetDofNumber(1)));
    }

    // ************************************************************************
    //      assemble and solve - TODO something like SolveStatic
    // ************************************************************************
    SimpleAssembler assembler(dofInfo);

    auto gradient = assembler.BuildVector(cellGroup, {displ}, GradientF);
    auto hessian = assembler.BuildMatrix(cellGroup, {displ}, Hessian0F);

    // build the constraint (modified) hessians and gradient
    // this should be moved the solver routine
    // solution = solver.solve(fullhessian, fullgradient, constraints);
    // no need to deal with modified vectors and matrices
    Eigen::SparseMatrix<double> hessianMod = cMatUnit.transpose() * hessian(displ, displ) * cMatUnit;
    Eigen::VectorXd residualMod = cMatUnit.transpose() * gradient[displ];
    Eigen::VectorXd deltaDisplacementsMod = EigenSparseSolve(hessianMod, residualMod, std::string("MumpsLDLT"));

    Eigen::VectorXd deltaDisplacements =
            cMatUnit * (-deltaDisplacementsMod) + constraints.GetSparseGlobalRhs(displ, numDofs, t);
    u[displ] = deltaDisplacements;

    // ************************************************************************
    //      merge dof values - TODO function MergeDofValues
    // ************************************************************************
    for (auto& node : mesh.NodesTotal(displ))
    {
        node.SetValue(0, u(displ, node.GetDofNumber(0)));
        node.SetValue(1, u(displ, node.GetDofNumber(1)));
    }

    Visualize::Visualizer visualize(cellGroup, Visualize::AverageHandler());
    visualize.DofValues(displ);

    auto stress = [linearElasticLaw, displ](const CellIpData& cellIpData) {
        return linearElasticLaw.Stress(cellIpData.Apply(displ, Nabla::Strain()));
    };
    visualize.CellData(stress, "Stress");

    visualize.CellData([](const CellIpData&) { return Eigen::Matrix<double, 1, 1>(7.0); }, "Seven");
    visualize.CellData([](const CellIpData& cipd) { return EigenCompanion::ToEigen(cipd.Ids().cellId); }, "CellId");
    visualize.WriteVtuFile("PatchTestVoronoi.vtu");

    Visualize::Visualizer pointVisualizer(cellGroup, Visualize::PointHandler(integrationType));
    pointVisualizer.DofValues(displ);
    pointVisualizer.CellData(stress, "Stress");
    pointVisualizer.WriteVtuFile("PatchTestPoint.vtu");

    // ************************************************************************
    //             check solution
    // ************************************************************************
    auto analyticDisplacementField = [=](Eigen::Vector2d coord) { return Eigen::Vector2d(coord[0] * 0.1, 0); };

    for (auto& node : geoMesh.NodesTotal())
    {
        auto coord = node.GetCoordinates();
        auto& displNode = mesh.NodeAtCoordinate(coord, displ);

        auto analyticSolution = analyticDisplacementField(coord);

        BOOST_TEST_MESSAGE("Node at " << coord.transpose() << " with dofs " << displNode.GetDofNumber(0) << ","
                                      << displNode.GetDofNumber(1));
        BOOST_TEST_MESSAGE("( " << displNode.GetValues().transpose() << " )");

        BoostUnitTest::CheckEigenMatrix(displNode.GetValues(), analyticSolution);
    }
}
