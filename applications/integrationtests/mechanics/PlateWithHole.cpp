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

Constraint::Constraints FixBottomAndLeft(MeshFem* rMesh, DofType disp)
{
    Constraint::Constraints constraints;

    auto nodesLeft = rMesh->NodesAtAxis(eDirection::X, disp, 0.);
    auto nodesBottom = rMesh->NodesAtAxis(eDirection::Y, disp, 0.);

    constraints.Add(disp, Constraint::Component(nodesLeft, {eDirection::X}));
    constraints.Add(disp, Constraint::Component(nodesBottom, {eDirection::Y}));

    return constraints;
}

BOOST_AUTO_TEST_CASE(PlateWithHole)
{
    auto binary = boost::unit_test::framework::master_test_suite().argv[0];
    boost::filesystem::path binaryPath = std::string(binary);
    binaryPath.remove_filename();

    std::string meshFile = binaryPath.string() + "/meshes/PlateWithHole.msh";

    auto meshGmsh = MeshGmsh(meshFile);
    auto& mesh = meshGmsh.GetMeshFEM();
    DofType disp("displacements", 2);
    AddDofInterpolation(&mesh, disp);

    const double E = 6174.;
    const double nu = 0.1415;
    Laws::LinearElastic<2> law(E, nu, ePlaneState::PLANE_STRESS);
    Integrands::MomentumBalance<2> momentumBalance(disp, law);

    Integrands::NeumannBc<2> neumannRight(disp, Test::PlateWithHoleAnalytical::PressureRight);
    Integrands::NeumannBc<2> neumannTop(disp, Test::PlateWithHoleAnalytical::PressureTop);

    auto Hessian0Plate = Bind(momentumBalance, &Integrands::MomentumBalance<2>::Hessian0);
    auto GradientPlate = Bind(momentumBalance, &Integrands::MomentumBalance<2>::Gradient);
    auto GradientRight = Bind(neumannRight, &Integrands::NeumannBc<2>::ExternalLoad);
    auto GradientTop = Bind(neumannTop, &Integrands::NeumannBc<2>::ExternalLoad);

    boost::ptr_vector<CellInterface> cells;
    Group<CellInterface> cellsPlate;
    IntegrationTypeTensorProduct<2> integrationTypeCells(2, eIntegrationMethod::GAUSS);
    for (auto& element : meshGmsh.GetPhysicalGroup("Plate"))
    {
        cells.push_back(new Cell(element, integrationTypeCells, 0));
        cellsPlate.Add(cells.back());
    }

    IntegrationTypeTensorProduct<1> integrationTypeBoundary(2, eIntegrationMethod::GAUSS);
    Group<CellInterface> cellsRight;
    for (auto& element : meshGmsh.GetPhysicalGroup("Right"))
    {
        cells.push_back(new Cell(element, integrationTypeBoundary, 0));
        cellsRight.Add(cells.back());
    }

    Group<CellInterface> cellsTop;
    for (auto& element : meshGmsh.GetPhysicalGroup("Top"))
    {
        cells.push_back(new Cell(element, integrationTypeBoundary, 0));
        cellsTop.Add(cells.back());
    }

    auto dispNodes = mesh.NodesTotal(disp);

    auto constraints = FixBottomAndLeft(&mesh, disp);
    auto dofInfo = DofNumbering::Build(dispNodes, disp, constraints);
    auto assembler = SimpleAssembler(dofInfo.numIndependentDofs, dofInfo.numDependentDofs);

    auto hessian0 = assembler.BuildMatrix(cellsPlate, {disp}, Hessian0Plate);
    auto gradient = assembler.BuildVector(cellsPlate, {disp}, GradientPlate);
    gradient += assembler.BuildVector(cellsRight, {disp}, GradientRight);
    gradient += assembler.BuildVector(cellsTop, {disp}, GradientTop);

    // std::cout << gradient.J[disp] << std::endl;
    std::cout << gradient.J[disp].norm() << std::endl;

    Eigen::VectorXd displ = -Eigen::SimplicialLLT<Eigen::SparseMatrix<double>>()
                                     .compute(hessian0.JJ(disp, disp))
                                     .solve(gradient.J[disp]);


    for (auto& node : dispNodes)
    {
        auto dofX = node.GetDofNumber(0);
        auto dofY = node.GetDofNumber(1);

        if (dofX < dofInfo.numIndependentDofs[disp])
            node.SetValue(0, displ[dofX]);
        if (dofY < dofInfo.numIndependentDofs[disp])
            node.SetValue(1, displ[dofY]);
    }
    gradient = assembler.BuildVector(cellsPlate, {disp}, GradientPlate);
    gradient += assembler.BuildVector(cellsRight, {disp}, GradientRight);
    gradient += assembler.BuildVector(cellsTop, {disp}, GradientTop);

    std::cout << gradient.J[disp].norm() << std::endl;

    Visualize::Visualizer<Visualize::VoronoiHandler> visu(cellsPlate, Visualize::VoronoiGeometryQuad(3));
    visu.DofValues(disp);
    auto analyticDisplacement = [E, nu](Eigen::Vector2d coords) {
        return Test::PlateWithHoleAnalytical::AnalyticDisplacement(coords, E, nu);
    };
    visu.PointData(analyticDisplacement, "analyticDisplacement");
    visu.WriteVtuFile("PlateWithHole.vtu");

    for (auto& element : meshGmsh.GetPhysicalGroup("Plate"))
    {
        const auto& coordElement = element.CoordinateElement();
        const auto& displElement = element.DofElement(disp);

        // we have isoparametric elements, so the nth coordinate node corresponds to the nth displacement node
        for (int iNode = 0; iNode < displElement.GetNumNodes(); ++iNode)
        {
            Eigen::Vector2d coords = coordElement.GetNode(iNode).GetValues();
            Eigen::Vector2d displs = displElement.GetNode(iNode).GetValues();
            Eigen::Vector2d correctDispl = Test::PlateWithHoleAnalytical::AnalyticDisplacement(coords, E, nu);

            BoostUnitTest::CheckEigenMatrix(displs, correctDispl, 1.e-5);
        }
    }
}
