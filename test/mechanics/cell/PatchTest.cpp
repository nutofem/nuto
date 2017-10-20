#include "BoostUnitTest.h"

#include <eigen3/Eigen/SparseLU>

#include "base/Group.h"

#include "mechanics/mesh/MeshFem.h"
#include "mechanics/mesh/MeshFemDofConvert.h"

#include "mechanics/interpolation/InterpolationQuadLinear.h"

#include "mechanics/integrationtypes/IntegrationTypeTensorProduct.h"

#include "mechanics/constitutive/laws/LinearElastic.h"

#include "mechanics/integrands/MomentumBalance.h"

#include "mechanics/cell/Cell.h"
#include "mechanics/cell/SimpleAssember.h"

using namespace NuTo;
using namespace NuTo::Groups;


MeshFem QuadPatchTestMesh()
{
    /*
     * Something like this:
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
    mesh.Elements.Add({{{n4, n4, n6, n7}, interpolation}});
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
    Group<NodeSimple> nodesConstrainedInX = Group<NodeSimple>(rMesh->NodeAtCoordinate(Eigen::Vector2d(0,0), dof));
    Group<NodeSimple> nodesConstrainedInY = rMesh->NodesAtAxis(eDirection::X, dof);


    Group<NodeSimple> nodesUnconstrainedInX = Difference(allNodes, nodesConstrainedInX);
    Group<NodeSimple> nodesUnconstrainedInY = Difference(allNodes, nodesConstrainedInY);

    int unconstrainedDofNumber = 0;
    int constrainedDofNumber = nodesUnconstrainedInX.Size() + nodesUnconstrainedInY.Size();

    for (auto& node : nodesUnconstrainedInX)
        node.SetDofNumber(0, unconstrainedDofNumber++);
    for (auto& node : nodesUnconstrainedInY)
        node.SetDofNumber(1, unconstrainedDofNumber++);
    for (auto& node : nodesConstrainedInX)
        node.SetDofNumber(0, constrainedDofNumber++);
    for (auto& node : nodesConstrainedInY)
        node.SetDofNumber(1, constrainedDofNumber++);


    DofInfo dofInfo;
    dofInfo.numDependentDofs[dof] = constrainedDofNumber;
    dofInfo.numIndependentDofs[dof] = unconstrainedDofNumber;
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
    //                 add boundary cells (TODO)
    // ************************************************************************
    Eigen::Vector2d pressureBC(42, 0);
    auto analyticDisplacementField = [=] (Eigen::Vector2d coord)
    {
        // deltaU / L = strain = E * stress = pressure/A
        const double strainX = E * pressureBC[0] / 10;
        const double strainY = E * nu * pressureBC[1] / 10;
        return Eigen::Vector2d(strainX * coord[0] / 10, strainY * coord[1] / 10 );
    };

    // ************************************************************************
    //                  assemble and solve
    // ************************************************************************
    SimpleAssembler assembler(dofInfo.numIndependentDofs, dofInfo.numDependentDofs);

    auto gradient = assembler.BuildVector(cellGroup, {&displ}, Integrands::TimeDependent::Gradient());
    auto hessian = assembler.BuildMatrix(cellGroup, {&displ}, Integrands::TimeDependent::Hessian0());

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver(hessian.JJ(displ, displ));
    Eigen::VectorXd newDisplacements = solver.solve(gradient.J[displ]);

    // merge
    int numUnconstrainedDofs = dofInfo.numIndependentDofs[displ];
    for(auto& node : mesh.NodesTotal(displ))
    {
        int dofX = node.GetDofNumber(0);
        int dofY = node.GetDofNumber(1);

        if (dofX < numUnconstrainedDofs)
            node.SetValue(0, newDisplacements[dofX]);
        if (dofY < numUnconstrainedDofs)
            node.SetValue(1, newDisplacements[dofY]);
    }


    for (auto& node : mesh.NodesTotal())
    {
        auto coord = node.GetValues();
        auto& displNode = mesh.NodeAtCoordinate(coord);

        auto analyticSolution = analyticDisplacementField(coord);
        BoostUnitTest::CheckEigenMatrix(displNode.GetValues(), analyticSolution);
    }
}

