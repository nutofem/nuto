#include "BoostUnitTest.h"
#include <boost/filesystem.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <string>
#include <iostream>

#include "base/Group.h"

#include "mechanics/dofs/DofNumbering.h"

#include "mechanics/mesh/MeshFem.h"
#include "mechanics/mesh/MeshGmsh.h"
#include "mechanics/mesh/MeshFemDofConvert.h"

#include "mechanics/interpolation/InterpolationQuadLinear.h"
#include "mechanics/interpolation/InterpolationTrussLinear.h"
#include "mechanics/integrationtypes/IntegrationTypeTensorProduct.h"

#include "mechanics/constraints/Constraints.h"
#include "mechanics/constraints/ConstraintCompanion.h"

#include "mechanics/constitutive/LinearElastic.h"

#include "mechanics/integrands/MomentumBalance.h"
#include "mechanics/integrands/NeumannBc.h"
#include "mechanics/integrands/Bind.h"

#include "mechanics/cell/Cell.h"
#include "mechanics/cell/SimpleAssembler.h"

#include "visualize/VoronoiHandler.h"
#include "visualize/VoronoiGeometries.h"
#include "visualize/Visualizer.h"

#include "PlateWithHoleAnalytic.h"

using namespace NuTo;

class CellAndAssemblyHelper
{
public:
    Group<CellInterface> AddCells(Group<ElementCollectionFem> elements, const IntegrationTypeBase& integrationType)
    {
        Group<CellInterface> cellGroup;
        for (auto& element : elements)
        {
            mCells.push_back(new Cell(element, integrationType, 0));
            cellGroup.Add(mCells.back());
        }
        return cellGroup;
    }

    void AddGradientFunction(Group<CellInterface> group, CellInterface::VectorFunction f)
    {
        mGradientFunctions.push_back({group, f});
    }

    void AddHessian0Function(Group<CellInterface> group, CellInterface::MatrixFunction f)
    {
        mHessian0Functions.push_back({group, f});
    }

    template <typename TSolver>
    void Solve(TSolver&& solver, Constraint::Constraints constraints, Group<NodeSimple> nodes, DofType dof)
    {
        auto dofInfo = DofNumbering::Build(nodes, dof, constraints);

        SimpleAssembler assembler(dofInfo.numIndependentDofs, dofInfo.numDependentDofs);

        auto hessianIterator = mHessian0Functions.begin();
        auto hessian0 = assembler.BuildMatrix(hessianIterator->first, {dof}, hessianIterator->second);
        hessianIterator++;
        for (; hessianIterator != mHessian0Functions.end(); ++hessianIterator)
            hessian0 += assembler.BuildMatrix(hessianIterator->first, {dof}, hessianIterator->second);

        auto gradientIterator = mGradientFunctions.begin();
        auto gradient = assembler.BuildVector(gradientIterator->first, {dof}, gradientIterator->second);
        gradientIterator++;
        for (; gradientIterator != mGradientFunctions.end(); ++gradientIterator)
            gradient += assembler.BuildVector(gradientIterator->first, {dof}, gradientIterator->second);

        Eigen::VectorXd newDisp = -solver.compute(hessian0.JJ(dof, dof)).solve(gradient.J[dof]);

        for (auto& node : nodes)
        {
            auto dofX = node.GetDofNumber(0);
            auto dofY = node.GetDofNumber(1);

            if (dofX < dofInfo.numIndependentDofs[dof])
                node.SetValue(0, newDisp[dofX]);
            if (dofY < dofInfo.numIndependentDofs[dof])
                node.SetValue(1, newDisp[dofY]);
        }

        gradientIterator = mGradientFunctions.begin();
        gradient = assembler.BuildVector(gradientIterator->first, {dof}, gradientIterator->second);
        gradientIterator++;
        for (; gradientIterator != mGradientFunctions.end(); ++gradientIterator)
            gradient += assembler.BuildVector(gradientIterator->first, {dof}, gradientIterator->second);

        double gradientNorm = gradient.J[dof].norm();
        if (gradientNorm > 1.e-10)
        {
            throw Exception(__PRETTY_FUNCTION__, "A single solve did not equilibrate the system.");
        }
    }


private:
    boost::ptr_vector<CellInterface> mCells;

    using GradientPair = std::pair<Group<CellInterface>, CellInterface::VectorFunction>;
    using Hessian0Pair = std::pair<Group<CellInterface>, CellInterface::MatrixFunction>;

    std::vector<GradientPair> mGradientFunctions;
    std::vector<Hessian0Pair> mHessian0Functions;
};


BOOST_AUTO_TEST_CASE(PlateWithHoleTest)
{
    auto binary = boost::unit_test::framework::master_test_suite().argv[0];
    boost::filesystem::path binaryPath = std::string(binary);
    binaryPath.remove_filename();

    std::string meshFile = binaryPath.string() + "/meshes/PlateWithHole.msh";
    MeshGmsh meshGmsh(meshFile);
    auto& mesh = meshGmsh.GetMeshFEM();
    DofType dof("Displacements", 2);
    AddDofInterpolation(&mesh, dof);

    CellAndAssemblyHelper helper;

    const double E = 6174.;
    const double nu = 0.1415;
    Laws::LinearElastic<2> law(E, nu, ePlaneState::PLANE_STRESS);
    Integrands::MomentumBalance<2> momentumBalance(dof, law);

    Integrands::NeumannBc<2> neumannRight(dof, Test::PlateWithHoleAnalytical::PressureRight);
    Integrands::NeumannBc<2> neumannTop(dof, Test::PlateWithHoleAnalytical::PressureTop);

    auto Hessian0Plate = Bind(momentumBalance, &Integrands::MomentumBalance<2>::Hessian0);
    auto GradientPlate = Bind(momentumBalance, &Integrands::MomentumBalance<2>::Gradient);
    auto GradientRight = Bind(neumannRight, &Integrands::NeumannBc<2>::ExternalLoad);
    auto GradientTop = Bind(neumannTop, &Integrands::NeumannBc<2>::ExternalLoad);

    IntegrationTypeTensorProduct<2> integrationTypeCells(2, eIntegrationMethod::GAUSS);
    IntegrationTypeTensorProduct<1> integrationTypeBoundary(2, eIntegrationMethod::GAUSS);
    auto cellsPlate = helper.AddCells(meshGmsh.GetPhysicalGroup("Plate"), integrationTypeCells);
    auto cellsRight = helper.AddCells(meshGmsh.GetPhysicalGroup("Right"), integrationTypeBoundary);
    auto cellsTop = helper.AddCells(meshGmsh.GetPhysicalGroup("Top"), integrationTypeBoundary);

    helper.AddHessian0Function(cellsPlate, Hessian0Plate);
    helper.AddGradientFunction(cellsPlate, GradientPlate);
    helper.AddGradientFunction(cellsRight, GradientRight);
    helper.AddGradientFunction(cellsTop, GradientTop);

    Constraint::Constraints constraints;
    auto nodesLeft = mesh.NodesAtAxis(eDirection::X, dof, 0.);
    auto nodesBottom = mesh.NodesAtAxis(eDirection::Y, dof, 0.);
    constraints.Add(dof, Constraint::Component(nodesLeft, {eDirection::X}));
    constraints.Add(dof, Constraint::Component(nodesBottom, {eDirection::Y}));

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    helper.Solve(solver, constraints, mesh.NodesTotal(dof), dof);


    Visualize::Visualizer<Visualize::VoronoiHandler> visu(cellsPlate, Visualize::VoronoiGeometryQuad(3));
    visu.DofValues(dof);
    auto analyticDisplacement = [E, nu](Eigen::Vector2d coords) {
        return Test::PlateWithHoleAnalytical::AnalyticDisplacement(coords, E, nu);
    };
    visu.PointData(analyticDisplacement, "analyticDisplacement");
    visu.WriteVtuFile("PlateWithHole.vtu");

    for (auto& element : meshGmsh.GetPhysicalGroup("Plate"))
    {
        const auto& coordElement = element.CoordinateElement();
        const auto& displElement = element.DofElement(dof);

        // we have isoparametric elements, so the nth coordinate node corresponds to the nth displacement node
        for (int iNode = 0; iNode < displElement.GetNumNodes(); ++iNode)
        {
            Eigen::Vector2d coords = coordElement.GetNode(iNode).GetValues();
            Eigen::Vector2d displs = displElement.GetNode(iNode).GetValues();
            Eigen::Vector2d correctDispl = analyticDisplacement(coords);

            BoostUnitTest::CheckEigenMatrix(displs, correctDispl, 1.e-5);
        }
    }
}
