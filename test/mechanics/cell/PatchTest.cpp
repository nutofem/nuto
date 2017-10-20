#include "BoostUnitTest.h"

#include <eigen3/Eigen/SparseLU>

#include "base/Group.h"

#include "mechanics/mesh/MeshFem.h"
#include "mechanics/mesh/MeshFemDofConvert.h"

#include "mechanics/interpolation/InterpolationQuadLinear.h"
#include "mechanics/interpolation/InterpolationTrussLinear.h"
#include "mechanics/integrationtypes/IntegrationTypeTensorProduct.h"

#include "mechanics/constitutive/laws/LinearElastic.h"

#include "mechanics/integrands/MomentumBalance.h"
#include "mechanics/integrands/NeumannBc.h"

#include "mechanics/cell/Cell.h"
#include "mechanics/cell/SimpleAssember.h"

using namespace NuTo;
using namespace NuTo::Groups;

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
    auto& n0 = mesh.Nodes.Add(Eigen::Vector2d(0, 0));
    auto& n1 = mesh.Nodes.Add(Eigen::Vector2d(10, 0));
    auto& n2 = mesh.Nodes.Add(Eigen::Vector2d(10, 10));
    auto& n3 = mesh.Nodes.Add(Eigen::Vector2d(0, 10));

    auto& n4 = mesh.Nodes.Add(Eigen::Vector2d(2, 2));
    auto& n5 = mesh.Nodes.Add(Eigen::Vector2d(8, 3));
    auto& n6 = mesh.Nodes.Add(Eigen::Vector2d(8, 7));
    auto& n7 = mesh.Nodes.Add(Eigen::Vector2d(4, 7));

    const auto& interpolation = mesh.CreateInterpolation(InterpolationQuadLinear(2));

    mesh.Elements.Add({{{n0, n1, n5, n4}, interpolation}});
    mesh.Elements.Add({{{n1, n2, n6, n5}, interpolation}});
    mesh.Elements.Add({{{n7, n6, n2, n3}, interpolation}});
    mesh.Elements.Add({{{n4, n5, n6, n7}, interpolation}});
    mesh.Elements.Add({{{n0, n4, n7, n3}, interpolation}});

    return mesh;
}

struct DofInfo
{
    DofContainer<int> numIndependentDofs;
    DofContainer<int> numDependentDofs;
};

DofInfo ManualDofNumbering(MeshFem* rMesh, DofType dof)
{
    // some manual dof numbering ...
    Group<NodeSimple> allNodes = rMesh->NodesTotal(dof);
    Group<NodeSimple> nodesConstrainedInX = rMesh->NodesAtAxis(eDirection::X, dof);
    Group<NodeSimple> nodesConstrainedInY = Group<NodeSimple>(rMesh->NodeAtCoordinate(Eigen::Vector2d(0, 0), dof));

    Group<NodeSimple> nodesUnconstrainedInX = Difference(allNodes, nodesConstrainedInX);
    Group<NodeSimple> nodesUnconstrainedInY = Difference(allNodes, nodesConstrainedInY);

    int dofNumber = 0;

    for (auto& node : nodesUnconstrainedInX)
        node.SetDofNumber(0, dofNumber++);
    for (auto& node : nodesUnconstrainedInY)
        node.SetDofNumber(1, dofNumber++);

    for (auto& node : nodesConstrainedInX)
        node.SetDofNumber(0, dofNumber++);
    for (auto& node : nodesConstrainedInY)
        node.SetDofNumber(1, dofNumber++);

    DofInfo dofInfo;
    dofInfo.numIndependentDofs[dof] = nodesUnconstrainedInX.Size() + nodesUnconstrainedInY.Size();
    dofInfo.numDependentDofs[dof] = nodesConstrainedInX.Size() + nodesConstrainedInY.Size();
    return dofInfo;
}

BOOST_AUTO_TEST_CASE(PatchTest)
{
    MeshFem mesh = QuadPatchTestMesh();
    DofType displ("displacements", 2);
    const auto& interpolation = mesh.CreateInterpolation(InterpolationQuadLinear(2));

    AddDofInterpolation(&mesh, displ, interpolation);
    DofInfo dofInfo = ManualDofNumbering(&mesh, displ);


    // ************************************************************************
    //                 add continuum cells
    // ************************************************************************
    constexpr double E = 20000;
    constexpr double nu = 0.2;
    Laws::LinearElastic<2> linearElasticLaw(E, nu, ePlaneState::PLANE_STRESS);
    Integrands::TimeDependent::MomentumBalance<2> momentumBalance(displ, linearElasticLaw);

    boost::ptr_vector<CellInterface> cellContainer;
    IntegrationTypeTensorProduct<2> integrationType(2, eIntegrationMethod::GAUSS);

    Group<CellInterface> cellGroup;
    for (auto& element : mesh.Elements)
    {
        cellContainer.push_back(new Cell(element, integrationType, momentumBalance));
        cellGroup.Add(cellContainer.back());
    }

    // ************************************************************************
    //                 add boundary cells
    // ************************************************************************

    // manually add the boundary element
    const auto& interpolationBc = mesh.CreateInterpolation(InterpolationTrussLinear(2));

    // extract existing nodes
    auto boundaryCoordNodes = mesh.NodesAtAxis(eDirection::X, 10);
    NodeSimple& nc1 = *boundaryCoordNodes.begin();
    NodeSimple& nc2 = *(boundaryCoordNodes.begin() + 1);

    auto boundaryDisplNodes = mesh.NodesAtAxis(eDirection::X, displ, 10);
    NodeSimple& nd1 = *boundaryDisplNodes.begin();
    NodeSimple& nd2 = *(boundaryDisplNodes.begin() + 1);

    // add the boundary element
    auto& boundaryElement = mesh.Elements.Add({{{nc1, nc2}, interpolationBc}});
    boundaryElement.AddDofElement(displ, {{nd1, nd2}, interpolationBc});

    IntegrationTypeTensorProduct<1> integrationTypeBc(1, eIntegrationMethod::GAUSS);
    Eigen::Vector2d pressureBC(1, 0);
    Integrands::TimeDependent::NeumannBc<2> neumannBc(displ, pressureBC);

    cellContainer.push_back(new Cell(boundaryElement, integrationTypeBc, neumannBc));
    cellGroup.Add(cellContainer.back());

    // ************************************************************************
    //                  assemble and solve
    // ************************************************************************
    SimpleAssembler assembler(dofInfo.numIndependentDofs, dofInfo.numDependentDofs);

    auto gradient = assembler.BuildVector(cellGroup, {&displ}, Integrands::TimeDependent::Gradient());
    auto hessian = assembler.BuildMatrix(cellGroup, {&displ}, Integrands::TimeDependent::Hessian0());

    Eigen::MatrixXd hessianDense(hessian.JJ(displ, displ));
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver(hessian.JJ(displ, displ));
    Eigen::VectorXd newDisplacements = solver.solve(gradient.J[displ]);

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
