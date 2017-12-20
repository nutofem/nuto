#include "BoostUnitTest.h"

#include <eigen3/Eigen/Dense> // for solve
#include "boost/ptr_container/ptr_vector.hpp"

#include "base/Group.h"

#include "mechanics/dofs/DofNumbering.h"

#include "mechanics/mesh/MeshFem.h"
#include "mechanics/mesh/MeshFemDofConvert.h"

#include "mechanics/interpolation/InterpolationQuadLinear.h"
#include "mechanics/interpolation/InterpolationTrussLinear.h"
#include "mechanics/integrationtypes/IntegrationTypeTensorProduct.h"

#include "mechanics/constraints/Constraints.h"
#include "mechanics/constraints/ConstraintCompanion.h"

#include "mechanics/constitutive/LinearElastic.h"

#include "mechanics/integrands/MomentumBalance.h"
#include "mechanics/integrands/NeumannBc.h"

#include "mechanics/cell/Cell.h"
#include "mechanics/cell/SimpleAssember.h"

using namespace NuTo;

//! @brief automatically create the lambda
//! [&](cellData, cellIpData) {return integrand.Gradient(cellData, cellIpData, 0); }
template <typename TObject, typename TReturn>
auto Bind(TObject& object, TReturn (TObject::*f)(const NuTo::CellData&, const NuTo::CellIpData&, double))
{
    return [&object, f](const NuTo::CellData& cellData, const NuTo::CellIpData& cellIpData) {
        return (object.*f)(cellData, cellIpData, /* deltaT = */ 0.);
    };
}
//! @brief automatically create the lambda
//! [&](cellData, cellIpData) {return integrand.Gradient(cellData, cellIpData); }
template <typename TObject, typename TReturn>
auto Bind(TObject& object, TReturn (TObject::*f)(const NuTo::CellData&, const NuTo::CellIpData&))
{
    return [&object, f](const NuTo::CellData& cellData, const NuTo::CellIpData& cellIpData) {
        return (object.*f)(cellData, cellIpData);
    };
}


MeshFem QuadPatchTestMesh()
{
    /* Something like this:
     *
     *    3-----------------------2
     * /| | - _        e2       / | -->
     * /| |     -7------------6   | -->
     * /| | e5  /     e4      |   | -->
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
    MeshFem mesh;
    NodeSimple& n0 = mesh.Nodes.Add(Eigen::Vector2d(0, 0));
    NodeSimple& n1 = mesh.Nodes.Add(Eigen::Vector2d(10, 0));
    NodeSimple& n2 = mesh.Nodes.Add(Eigen::Vector2d(10, 10));
    NodeSimple& n3 = mesh.Nodes.Add(Eigen::Vector2d(0, 10));

    NodeSimple& n4 = mesh.Nodes.Add(Eigen::Vector2d(2, 2));
    NodeSimple& n5 = mesh.Nodes.Add(Eigen::Vector2d(8, 3));
    NodeSimple& n6 = mesh.Nodes.Add(Eigen::Vector2d(8, 7));
    NodeSimple& n7 = mesh.Nodes.Add(Eigen::Vector2d(4, 7));

    const InterpolationSimple& interpolation = mesh.CreateInterpolation(InterpolationQuadLinear());

    mesh.Elements.Add({{{n0, n1, n5, n4}, interpolation}});
    mesh.Elements.Add({{{n1, n2, n6, n5}, interpolation}});
    mesh.Elements.Add({{{n7, n6, n2, n3}, interpolation}});
    mesh.Elements.Add({{{n4, n5, n6, n7}, interpolation}});
    mesh.Elements.Add({{{n0, n4, n7, n3}, interpolation}});

    return mesh;
}

Constraint::Constraints DefineConstraints(MeshFem* rMesh, DofType dof)
{
    Constraint::Constraints constraints;

    Group<NodeSimple> nodesConstrainedInX = rMesh->NodesAtAxis(eDirection::X, dof);
    Group<NodeSimple> nodesConstrainedInY = Group<NodeSimple>(rMesh->NodeAtCoordinate(Eigen::Vector2d(0, 0), dof));

    constraints.Add(dof, Constraint::Component(nodesConstrainedInX, {eDirection::X}));
    constraints.Add(dof, Constraint::Component(nodesConstrainedInY, {eDirection::Y}));

    return constraints;
}

BOOST_AUTO_TEST_CASE(PatchTestForce)
{
    MeshFem mesh = QuadPatchTestMesh();
    DofType displ("displacements", 2);
    const InterpolationSimple& interpolation = mesh.CreateInterpolation(InterpolationQuadLinear());

    AddDofInterpolation(&mesh, displ, interpolation);

    Constraint::Constraints constraints = DefineConstraints(&mesh, displ);
    DofNumbering::DofInfo dofInfo = DofNumbering::Build(mesh.NodesTotal(displ), displ, constraints);


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
    for (ElementCollection& element : mesh.Elements)
    {
        cellContainer.push_back(new Cell(element, integrationType, cellId++));
        momentumBalanceCells.Add(cellContainer.back());
    }

    // ************************************************************************
    //                 add boundary cells
    // ************************************************************************

    // manually add the boundary element
    const InterpolationSimple& interpolationBc = mesh.CreateInterpolation(InterpolationTrussLinear());

    // extract existing nodes
    Group<NodeSimple> boundaryCoordNodes = mesh.NodesAtAxis(eDirection::X, 10);
    NodeSimple& nc1 = *boundaryCoordNodes.begin();
    NodeSimple& nc2 = *(boundaryCoordNodes.begin() + 1);

    Group<NodeSimple> boundaryDisplNodes = mesh.NodesAtAxis(eDirection::X, displ, 10);
    NodeSimple& nd1 = *boundaryDisplNodes.begin();
    NodeSimple& nd2 = *(boundaryDisplNodes.begin() + 1);

    // add the boundary element
    ElementCollectionFem& boundaryElement = mesh.Elements.Add({{{nc1, nc2}, interpolationBc}});
    boundaryElement.AddDofElement(displ, {{nd1, nd2}, interpolationBc});

    IntegrationTypeTensorProduct<1> integrationTypeBc(1, eIntegrationMethod::GAUSS);
    Eigen::Vector2d pressureBC(1, 0);
    Integrands::NeumannBc<2> neumannBc(displ, pressureBC);
    auto NeumannLoad = Bind(neumannBc, &Integrands::NeumannBc<2>::ExternalLoad);

    cellContainer.push_back(new Cell(boundaryElement, integrationTypeBc, cellId++));
    auto& neumannCell = cellContainer.back();

    // ************************************************************************
    //                  assemble and solve
    // ************************************************************************
    SimpleAssembler assembler(dofInfo.numIndependentDofs, dofInfo.numDependentDofs);

    GlobalDofVector gradient = assembler.BuildVector(momentumBalanceCells, {displ}, MomentumGradientF);
    gradient += assembler.BuildVector({neumannCell}, {displ}, NeumannLoad);

    GlobalDofMatrixSparse hessian = assembler.BuildMatrix(momentumBalanceCells, {displ}, MomentumHessian0F);
    // no hessian for the neumann bc integrand (external load)

    Eigen::MatrixXd hessianDense(hessian.JJ(displ, displ));
    Eigen::VectorXd newDisplacements = hessianDense.ldlt().solve(gradient.J[displ]);

    // merge dof values
    int numUnconstrainedDofs = dofInfo.numIndependentDofs[displ];
    for (auto& node : mesh.NodesTotal(displ))
    {
        int dofX = node.GetDofNumber(0);
        int dofY = node.GetDofNumber(1);

        if (dofX < numUnconstrainedDofs)
            node.SetValue(0, newDisplacements[dofX]);
        if (dofY < numUnconstrainedDofs)
            node.SetValue(1, newDisplacements[dofY]);
    }

    auto analyticDisplacementField = [=](Eigen::Vector2d coord) {
        // pressure / A = sigma = E * strain = E * deltaU / L
        return Eigen::Vector2d(pressureBC[0] / E * coord[0], -nu * pressureBC[0] / E * coord[1]);
    };

    for (NodeSimple& node : mesh.NodesTotal())
    {
        Eigen::VectorXd coord = node.GetValues();
        NodeSimple& displNode = mesh.NodeAtCoordinate(coord, displ);

        Eigen::VectorXd analyticSolution = analyticDisplacementField(coord);

        BOOST_TEST_MESSAGE("Node at " << coord.transpose() << " with dofs " << displNode.GetDofNumber(0) << ","
                                      << displNode.GetDofNumber(1));
        BOOST_TEST_MESSAGE("( " << displNode.GetValues().transpose() << " )");

        BoostUnitTest::CheckEigenMatrix(displNode.GetValues(), analyticSolution);
    }
}

BOOST_AUTO_TEST_CASE(PatchTestDispl)
{
    MeshFem mesh = QuadPatchTestMesh();
    DofType displ("displacements", 2);
    const auto& interpolation = mesh.CreateInterpolation(InterpolationQuadLinear());

    AddDofInterpolation(&mesh, displ, interpolation);

    auto constraints = DefineConstraints(&mesh, displ); // fixed boundary conditions
    Group<NodeSimple> rightBoundary = mesh.NodesAtAxis(eDirection::X, displ, 10);
    const double boundaryDisplacement = 1.;
    constraints.Add(displ, Constraint::Component(rightBoundary, {eDirection::X}, boundaryDisplacement));

    DofNumbering::DofInfo dofInfo = DofNumbering::Build(mesh.NodesTotal(displ), displ, constraints);
    const int numDofs = dofInfo.numIndependentDofs[displ] + dofInfo.numDependentDofs[displ];
    const int numDepDofs = dofInfo.numDependentDofs[displ];
    Eigen::MatrixXd CMat = constraints.BuildConstraintMatrix(displ, numDofs - numDepDofs);


    BOOST_TEST_MESSAGE("CMat \n" << CMat);

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
    for (auto& element : mesh.Elements)
    {
        cellContainer.push_back(new Cell(element, integrationType, cellId++));
        cellGroup.Add(cellContainer.back());
    }

    // ************************************************************************
    //      assemble and solve - TODO something like SolveStatic
    // ************************************************************************
    SimpleAssembler assembler(dofInfo.numIndependentDofs, dofInfo.numDependentDofs);

    auto gradient = assembler.BuildVector(cellGroup, {displ}, GradientF);
    auto hessian = assembler.BuildMatrix(cellGroup, {displ}, Hessian0F);

    Eigen::MatrixXd kJJ = Eigen::MatrixXd(hessian.JJ(displ, displ));
    Eigen::MatrixXd kJK = Eigen::MatrixXd(hessian.JK(displ, displ));
    Eigen::MatrixXd kKJ = Eigen::MatrixXd(hessian.KJ(displ, displ));
    Eigen::MatrixXd kKK = Eigen::MatrixXd(hessian.KK(displ, displ));
    BOOST_TEST_MESSAGE("hessian JJ \n" << kJJ);
    BOOST_TEST_MESSAGE("hessian JK \n" << kJK);
    BOOST_TEST_MESSAGE("hessian KJ \n" << kKJ);
    BOOST_TEST_MESSAGE("hessian KK \n" << kKK);


    BOOST_TEST_MESSAGE("GradientJ \n" << gradient.J);
    BOOST_TEST_MESSAGE("GradientK \n" << gradient.K);

    //      have a look at DISS_UNGER, page 28 for all that CMat stuff.
    Eigen::MatrixXd Kmod = kJJ - CMat.transpose() * kKJ - kJK * CMat + CMat.transpose() * kKK * CMat;
    Eigen::VectorXd Rmod = gradient.J[displ] - CMat.transpose() * gradient.K[displ];
    Eigen::VectorXd RmodConstrained = (kJK - CMat.transpose() * kKK) * (-constraints.GetRhs(displ, 0));

    Eigen::VectorXd newDisplacementsJ = Kmod.ldlt().solve(Rmod + RmodConstrained);
    Eigen::VectorXd newDisplacementsK = -CMat * newDisplacementsJ + constraints.GetRhs(displ, 0);

    BOOST_TEST_MESSAGE("DeltaD J \n" << newDisplacementsJ);
    BOOST_TEST_MESSAGE("DeltaD K \n" << newDisplacementsK);

    // ************************************************************************
    //      merge dof values - TODO function MergeDofValues
    // ************************************************************************
    int numUnconstrainedDofs = dofInfo.numIndependentDofs[displ];
    for (auto& node : mesh.NodesTotal(displ))
    {
        int dofX = node.GetDofNumber(0);
        int dofY = node.GetDofNumber(1);

        if (dofX < numUnconstrainedDofs)
            node.SetValue(0, newDisplacementsJ[dofX]);
        else
            node.SetValue(0, newDisplacementsK[dofX - numUnconstrainedDofs]);


        if (dofY < numUnconstrainedDofs)
            node.SetValue(1, newDisplacementsJ[dofY]);
        else
            node.SetValue(1, newDisplacementsK[dofY - numUnconstrainedDofs]);
    }


    // ************************************************************************
    //             check solution
    // ************************************************************************
    auto analyticDisplacementField = [=](Eigen::Vector2d coord) { return Eigen::Vector2d(coord[0] * 0.1, 0); };

    for (auto& node : mesh.NodesTotal())
    {
        auto coord = node.GetValues();
        auto& displNode = mesh.NodeAtCoordinate(coord, displ);

        auto analyticSolution = analyticDisplacementField(coord);

        BOOST_TEST_MESSAGE("Node at " << coord.transpose() << " with dofs " << displNode.GetDofNumber(0) << ","
                                      << displNode.GetDofNumber(1));
        BOOST_TEST_MESSAGE("( " << displNode.GetValues().transpose() << " )");

        BoostUnitTest::CheckEigenMatrix(displNode.GetValues(), analyticSolution);
    }
}
