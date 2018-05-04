#include "BoostUnitTest.h"
#include <type_traits>
#include <iostream>

#include "nuto/mechanics/iga/Nurbs.h"
#include "nuto/mechanics/elements/ElementIga.h"
#include "nuto/mechanics/dofs/DofType.h"

#include "nuto/base/Group.h"
#include "nuto/base/ValueVector.h"
#include "nuto/mechanics/DirectionEnum.h"
#include "nuto/mechanics/nodes/NodeSimple.h"
#include "nuto/mechanics/elements/ElementCollection.h"

#include "nuto/mechanics/constraints/Constraints.h"
#include "nuto/mechanics/constraints/ConstraintCompanion.h"

#include "nuto/mechanics/dofs/DofNumbering.h"

#include "nuto/mechanics/constitutive/LinearElastic.h"

#include "nuto/mechanics/integrands/MomentumBalance.h"
#include "nuto/mechanics/integrands/NeumannBc.h"
#include "nuto/mechanics/integrands/Bind.h"

#include "nuto/mechanics/cell/Cell.h"
#include "nuto/mechanics/cell/SimpleAssembler.h"

#include "boost/ptr_container/ptr_vector.hpp"

#include "nuto/mechanics/integrationtypes/IntegrationTypeTensorProduct.h"

#include "nuto/mechanics/solver/Solve.h"

#include "nuto/visualize/AverageGeometries.h"
#include "nuto/visualize/VoronoiGeometries.h"
#include "nuto/visualize/PointHandler.h"
#include "nuto/visualize/Visualizer.h"
#include "nuto/visualize/VoronoiGeometries.h"
#include "nuto/math/EigenCompanion.h"

#include <Eigen/Dense> // for determinant

#include "nuto/mechanics/mesh/MeshIga.h"
using namespace NuTo;

void QuadPatchTestNurbs(MeshIga<2>& mesh, const NuTo::DofType& coord, const NuTo::DofType& displ, int& num)
{
    std::vector<NodeSimple*> controlPointsCoordinates;

    num = 3;
    double incr = 1. / (num - 1);
    for (int i = 0; i < num; i++)
        for (int j = 0; j < num; j++)
        {
            auto& node = mesh.mNodes.Add(Eigen::Vector2d({j * incr, i * incr}));
            controlPointsCoordinates.push_back(&node);
        }

    std::vector<NodeSimple*> controlPointsDispl;
    for (int i = 0; i < num; i++)
        for (int j = 0; j < num; j++)
        {
            auto& node = mesh.mNodes.Add(Eigen::Vector2d({0, 0}));
            controlPointsDispl.push_back(&node);
        }

    std::vector<double> knots1 = {0, 0, 0, 1, 1, 1};
    std::vector<double> knots2 = {0, 0, 0, 1, 1, 1};
    std::array<std::vector<double>, 2> knots = {{knots1, knots2}};

    std::vector<double> weights = {1, 1, 1, 1, 1, 1, 1, 1, 1};

    std::array<int, 2> degree = {{2, 2}};

    mesh.mDofInterpolations.Insert(coord, Nurbs<2>(knots, controlPointsCoordinates, weights, degree));

    mesh.mDofInterpolations.Insert(displ, Nurbs<2>(knots, controlPointsDispl, weights, degree));

    std::array<int, 2> knotIDsCell = {{2, 2}};
    ElementIga<2> igaCoordinates(knotIDsCell, mesh.mDofInterpolations.At(coord));
    ElementIga<2> igaDispl(knotIDsCell, mesh.mDofInterpolations.At(displ));

    mesh.mElements.Add(igaCoordinates);
    mesh.mElements[0].AddDofElement(displ, igaDispl);
}

void QuadPatchTestNurbsFourElements(MeshIga<2>& mesh, const NuTo::DofType& coord, const NuTo::DofType& displ, int& num)
{
    std::vector<NodeSimple*> controlPointsCoordinates;

    num = 4;
    double incr = 1. / (num - 1);
    for (int i = 0; i < num; i++)
        for (int j = 0; j < num; j++)
        {
            auto& node = mesh.mNodes.Add(Eigen::Vector2d({j * incr, i * incr}));
            controlPointsCoordinates.push_back(&node);
            BOOST_TEST_MESSAGE(Eigen::Vector2d({j * incr, i * incr}).transpose());
        }

    std::vector<NodeSimple*> controlPointsDispl;
    for (int i = 0; i < num; i++)
        for (int j = 0; j < num; j++)
        {
            auto& node = mesh.mNodes.Add(Eigen::Vector2d({0, 0}));
            controlPointsDispl.push_back(&node);
        }

    std::vector<double> knots1 = {0, 0, 0, 0.5, 1, 1, 1};
    std::vector<double> knots2 = {0, 0, 0, 0.5, 1, 1, 1};
    std::array<std::vector<double>, 2> knots = {{knots1, knots2}};

    std::vector<double> weights = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    std::array<int, 2> degree = {{2, 2}};

    Nurbs<2> nurbs_coord(knots, controlPointsCoordinates, weights, degree);

    std::vector<double> ips;
    ips.push_back(-0.577350269189626);
    ips.push_back(0.);
    ips.push_back(+0.577350269189626);
    for (double value : ips)
    {
        Eigen::VectorXd parameter(2);
        parameter << 1, value;
        //        Eigen::VectorXd coordinateEvaluation = nurbs_coord.Evaluate(parameter, 0);
        //        BOOST_TEST_MESSAGE("coordinateEvaluation \n" << coordinateEvaluation);

        std::array<int, 2> knotID = {{3, 2}};
        const std::array<int, 2> spanIdx = nurbs_coord.FindSpan(nurbs_coord.Transformation(parameter, knotID));
        Eigen::MatrixXd mat = nurbs_coord.BasisFunctionsAndDerivativesRational(
                1, nurbs_coord.Transformation(parameter, knotID), spanIdx);
        Eigen::VectorXd coords = nurbs_coord.GetControlPointCoordinatesElement(knotID);

        Eigen::MatrixXd nodeBlockCoordinates =
                Eigen::Map<Eigen::Matrix<double, 2, Eigen::Dynamic>>(coords.data(), 2, mat.rows());

        BOOST_TEST_MESSAGE("2D \n" << (nodeBlockCoordinates * mat).col(1).norm());
    }

    mesh.mDofInterpolations.Insert(coord, nurbs_coord);
    mesh.mDofInterpolations.Insert(displ, Nurbs<2>(knots, controlPointsDispl, weights, degree));

    for (int i = 2; i <= num - 1; i++)
        for (int j = 2; j <= num - 1; j++)
        {
            std::array<int, 2> knotIDsCell = {{i, j}};
            ElementIga<2> igaCoordinates(knotIDsCell, mesh.mDofInterpolations.At(coord));
            ElementIga<2> igaDispl(knotIDsCell, mesh.mDofInterpolations.At(displ));

            ElementCollectionIga<2>& elementcollection = mesh.mElements.Add(igaCoordinates);
            elementcollection.AddDofElement(displ, igaDispl);
        }
}

// BOOST_AUTO_TEST_CASE(PatchTestIgaDispl)
//{
//    MeshIga<2> mesh;
//    NuTo::DofType coord("coordinates", 2);
//    NuTo::DofType displ("displacements", 2);
//    int num = 0;
//    QuadPatchTestNurbsFourElements(mesh, coord, displ, num);
//    // QuadPatchTestNurbs(mesh, coord, displ, num);

//    // %%%%%%%%%%%%%%%%%%%% Constraints %%%%%%%%%%%%%%%%%%%%%% //
//    NuTo::Constraint::Constraints constraints;

//    NuTo::Group<NuTo::NodeSimple> nodesConstrained;
//    for (int i = 0; i < num; i++)
//    {
//        Eigen::Vector2i nodeIds(0, i);
//        nodesConstrained.Add(*mesh.mDofInterpolations.At(displ).GetControlPoint(nodeIds));
//    }
//    constraints.Add(displ, NuTo::Constraint::Component(nodesConstrained, {NuTo::eDirection::X}));

//    NuTo::Group<NuTo::NodeSimple> nodesConstrained1;
//    Eigen::Vector2i nodeIds(0, 0);
//    nodesConstrained1.Add(*mesh.mDofInterpolations.At(displ).GetControlPoint(nodeIds));
//    constraints.Add(displ, NuTo::Constraint::Component(nodesConstrained1, {NuTo::eDirection::Y}));

//    NuTo::Group<NuTo::NodeSimple> nodesDisplConstrained;
//    for (int i = 0; i < num; i++)
//    {
//        Eigen::Vector2i nodeIds(num - 1, i);
//        nodesDisplConstrained.Add(*mesh.mDofInterpolations.At(displ).GetControlPoint(nodeIds));
//    }
//    constraints.Add(displ, NuTo::Constraint::Component(nodesDisplConstrained, {NuTo::eDirection::X}, 1));

//    NuTo::Group<NuTo::NodeSimple> nodesDispl;
//    for (int i = 0; i < num; i++)
//        for (int j = 0; j < num; j++)
//        {
//            nodeIds << i, j;
//            nodesDispl.Add(*mesh.mDofInterpolations.At(displ).GetControlPoint(nodeIds));
//        }

//    NuTo::DofInfo dofInfo = NuTo::DofNumbering::Build(nodesDispl, displ, constraints);

//    // %%%%%%%%%%%%%%%%%%%%%%%% material and integrand %%%%%%%%%%%%%%%%%%%%%%%%%%%% //
//    constexpr double E = 20000;
//    constexpr double nu = 0.0;
//    NuTo::Laws::LinearElastic<2> linearElasticLaw(E, nu, NuTo::ePlaneState::PLANE_STRESS);
//    Integrands::MomentumBalance<2> momentumBalance(displ, linearElasticLaw);
//    auto GradientF = Bind(momentumBalance, &Integrands::MomentumBalance<2>::Gradient);
//    auto Hessian0F = Bind(momentumBalance, &Integrands::MomentumBalance<2>::Hessian0);

//    boost::ptr_vector<CellInterface> cellContainer;
//    int numIPs = 3;
//    IntegrationTypeTensorProduct<2> integrationType(numIPs, eIntegrationMethod::GAUSS);

//    Group<CellInterface> cellGroup;
//    int cellId = 0;
//    for (auto& element : mesh.mElements)
//    {
//        cellContainer.push_back(new Cell(element, integrationType, cellId++));
//        cellGroup.Add(cellContainer.back());
//    }

//    // %%%%%%%%%%%%%%%%%%%%%%%% assemble and solve %%%%%%%%%%%%%%%%%%%%%%%%%%%% //
//    SimpleAssembler assembler(dofInfo);

//    auto gradient = assembler.BuildVector(cellGroup, {displ}, GradientF);
//    auto hessian = assembler.BuildMatrix(cellGroup, {displ}, Hessian0F);

//    int numIndependentDofs = dofInfo.numIndependentDofs[displ];

//    Eigen::MatrixXd CMat = constraints.BuildConstraintMatrix(displ, numIndependentDofs);

//    BOOST_TEST_MESSAGE("CMat \n" << CMat);

//    GlobalDofVector u = Solve(hessian, gradient, constraints, displ, numIndependentDofs, 0.0);

//    // %%%%%%%%%%%%%%%%% merge dof values %%%%%%%%% //
//    for (auto& node : nodesDispl)
//    {
//        node.SetValue(0, u(displ, node.GetDofNumber(0)));
//        node.SetValue(1, u(displ, node.GetDofNumber(1)));
//    }

//    //    Visualize::Spacing s = Visualize::Spacing::GAUSS;
//    Visualize::Visualizer visualize(cellGroup, Visualize::VoronoiHandler(Visualize::VoronoiGeometryQuad(numIPs)));
//    visualize.DofValues(displ);

//    auto stress = [linearElasticLaw, displ](const CellIpData& cellIpData) {
//        return linearElasticLaw.Stress(cellIpData.Apply(displ, Nabla::Strain()));
//    };
//    visualize.CellData(stress, "Stress");

//    //    visualize.CellData([](const CellIpData&) { return EigenCompanion::ToEigen(7.0); }, "Seven");
//    visualize.CellData([](const CellIpData& cipd) { return EigenCompanion::ToEigen(cipd.Ids().cellId); }, "CellId");
//    visualize.WriteVtuFile("voronoiIgaDirichlet.vtu");
//}


// BOOST_AUTO_TEST_CASE(PatchTestIgaNeumann)
//{
//    MeshIga<2> mesh;
//    NuTo::DofType coord("coordinates", 2);
//    NuTo::DofType displ("displacements", 2);
//    int num = 0;
//    QuadPatchTestNurbsFourElements(mesh, coord, displ, num);
//    // QuadPatchTestNurbs(mesh, coord, displ, num);

//    // %%%%%%%%%%%%%%%%%%%% Constraints %%%%%%%%%%%%%%%%%%%%%% //

//    NuTo::Constraint::Constraints constraints;

//    NuTo::Group<NuTo::NodeSimple> nodesConstrained;
//    for (int i = 0; i < num; i++)
//    {
//        Eigen::Vector2i nodeIds(0, i);
//        nodesConstrained.Add(*mesh.mDofInterpolations.At(displ).GetControlPoint(nodeIds));
//    }
//    constraints.Add(displ, NuTo::Constraint::Component(nodesConstrained, {NuTo::eDirection::X}));

//    NuTo::Group<NuTo::NodeSimple> nodesConstrained1;
//    Eigen::Vector2i nodeIds(0, 0);
//    nodesConstrained1.Add(*mesh.mDofInterpolations.At(displ).GetControlPoint(nodeIds));
//    constraints.Add(displ, NuTo::Constraint::Component(nodesConstrained1, {NuTo::eDirection::Y}));

//    // %%%%%%%%%%%%%%%%%%%% All nodes aka. CPs %%%%%%%%%%%%%%%%%%%%%% //

//    NuTo::Group<NuTo::NodeSimple> nodesDispl;
//    for (int i = 0; i < num; i++)
//        for (int j = 0; j < num; j++)
//        {
//            nodeIds << i, j;
//            nodesDispl.Add(*mesh.mDofInterpolations.At(displ).GetControlPoint(nodeIds));
//        }

//    NuTo::DofInfo dofInfo = NuTo::DofNumbering::Build(nodesDispl, displ, constraints);

//    // %%%%%%%%%%%%%%%%%%%%%%%% material and integrand %%%%%%%%%%%%%%%%%%%%%%%%%%%% //
//    constexpr double E = 20000;
//    constexpr double nu = 0.0;
//    NuTo::Laws::LinearElastic<2> linearElasticLaw(E, nu, NuTo::ePlaneState::PLANE_STRESS);
//    Integrands::MomentumBalance<2> momentumBalance(displ, linearElasticLaw);
//    auto GradientF = Bind(momentumBalance, &Integrands::MomentumBalance<2>::Gradient);
//    auto Hessian0F = Bind(momentumBalance, &Integrands::MomentumBalance<2>::Hessian0);

//    boost::ptr_vector<CellInterface> cellContainer;
//    int numIPs = 3;
//    IntegrationTypeTensorProduct<2> integrationType(numIPs, eIntegrationMethod::GAUSS);

//    Group<CellInterface> cellGroup;
//    int cellId = 0;
//    for (auto& element : mesh.mElements)
//    {
//        cellContainer.push_back(new Cell(element, integrationType, cellId++));
//        cellGroup.Add(cellContainer.back());
//    }

//    // %%%%%%%%%%%%%%%%%%%% Neumann BC %%%%%%%%%%%%%%%%%%%%%% //
//    std::vector<NodeSimple*> controlPointsCoordinatesNBC;
//    std::vector<NodeSimple*> controlPointsDisplNBC;

//    for (int i = 0; i < num; i++)
//    {
//        Eigen::Vector2i nodeIds(num - 1, i);
//        controlPointsCoordinatesNBC.push_back(mesh.mDofInterpolations.At(coord).GetControlPoint(nodeIds));
//        controlPointsDisplNBC.push_back(mesh.mDofInterpolations.At(displ).GetControlPoint(nodeIds));
//    }

//    std::vector<double> knots1 = {0, 0, 0, 0.5, 1, 1, 1};
//    std::vector<double> weightsNBC = {1, 1, 1, 1};
//    std::array<int, 1> degreeNBC = {{2}};

//    Nurbs<1> nurbsNBC_coord({{knots1}}, controlPointsCoordinatesNBC, weightsNBC, degreeNBC);

//    std::vector<double> ips;
//    ips.push_back(-0.577350269189626);
//    ips.push_back(0.);
//    ips.push_back(+0.577350269189626);
//    for (double value : ips)
//    {
//        Eigen::VectorXd parameter(1);
//        parameter << value;
//        //        Eigen::VectorXd coordinateEvaluation = nurbsNBC_coord.Evaluate(parameter, 0);
//        //        BOOST_TEST_MESSAGE("coordinateEvaluation \n" << coordinateEvaluation);

//        std::array<int, 1> knotID = {{2}};
//        const std::array<int, 1> spanIdx = nurbsNBC_coord.FindSpan(nurbsNBC_coord.Transformation(parameter, knotID));

//        Eigen::VectorXd coords = nurbsNBC_coord.GetControlPointCoordinatesElement(knotID);
//        Eigen::MatrixXd mat = nurbsNBC_coord.BasisFunctionsAndDerivativesRational(
//                1, nurbsNBC_coord.Transformation(parameter, knotID), spanIdx);

//        Eigen::MatrixXd nodeBlockCoordinates =
//                Eigen::Map<Eigen::Matrix<double, 2, Eigen::Dynamic>>(coords.data(), 2, mat.rows());

//        BOOST_TEST_MESSAGE("1D \n" << (nodeBlockCoordinates * mat).norm());
//    }

//    Nurbs<1> nurbsNBC_disp({{knots1}}, controlPointsDisplNBC, weightsNBC, degreeNBC);

//    std::vector<ElementCollectionIga<1>> vectorElementsNBC;
//    for (int i = 2; i <= num - 1; i++)
//    {
//        std::array<int, 1> knotIDsCell = {{i}};
//        ElementIga<1> igaCoordinates(knotIDsCell, nurbsNBC_coord);
//        ElementIga<1> igaDispl(knotIDsCell, nurbsNBC_disp);

//        vectorElementsNBC.push_back(igaCoordinates);
//        vectorElementsNBC.back().AddDofElement(displ, igaDispl);
//    }

//    IntegrationTypeTensorProduct<1> integrationTypeBc(3, eIntegrationMethod::GAUSS);
//    Eigen::Vector2d pressureBC(10, 0);
//    Integrands::NeumannBc<2> neumannBc(displ, pressureBC);
//    auto NeumannLoad = Bind(neumannBc, &Integrands::NeumannBc<2>::ExternalLoad);

//    NuTo::Group<CellInterface> neumannCells;

//    cellContainer.push_back(new Cell(vectorElementsNBC[0], integrationTypeBc, cellId++));
//    neumannCells.Add(cellContainer.back());
//    cellContainer.push_back(new Cell(vectorElementsNBC[1], integrationTypeBc, cellId++));
//    neumannCells.Add(cellContainer.back());

//    // %%%%%%%%%%%%%%%%%%%%%%%% assemble and solve %%%%%%%%%%%%%%%%%%%%%%%%%%%% //
//    SimpleAssembler assembler(dofInfo);

//    auto gradient = assembler.BuildVector(neumannCells, {displ}, NeumannLoad);
//    gradient += assembler.BuildVector(cellGroup, {displ}, GradientF);

//    auto hessian = assembler.BuildMatrix(cellGroup, {displ}, Hessian0F);

//    int numIndependentDofs = dofInfo.numIndependentDofs[displ];

//    Eigen::MatrixXd CMat = constraints.BuildConstraintMatrix(displ, numIndependentDofs);

//    BOOST_TEST_MESSAGE("CMat \n" << CMat);

//    GlobalDofVector u = Solve(hessian, gradient, constraints, displ, numIndependentDofs, 0.0);

//    // %%%%%%%%%%%%%%%%% merge dof values %%%%%%%%% //
//    for (auto& node : nodesDispl)
//    {
//        node.SetValue(0, u(displ, node.GetDofNumber(0)));
//        node.SetValue(1, u(displ, node.GetDofNumber(1)));
//    }

//    //    Visualize::Spacing s = Visualize::Spacing::GAUSS;
//    Visualize::Visualizer visualize(cellGroup, Visualize::VoronoiHandler(Visualize::VoronoiGeometryQuad(numIPs)));
//    visualize.DofValues(displ);

//    auto stress = [linearElasticLaw, displ](const CellIpData& cellIpData) {
//        return linearElasticLaw.Stress(cellIpData.Apply(displ, Nabla::Strain()));
//    };
//    visualize.CellData(stress, "Stress");

//    visualize.CellData([](const CellIpData& cipd) { return EigenCompanion::ToEigen(cipd.Ids().cellId); }, "CellId");
//    visualize.WriteVtuFile("voronoiIgaNeumann.vtu");
//}

void PlateWithHoleGeometry(MeshIga<2>& mesh, const NuTo::DofType& coord, const NuTo::DofType& displ, int refinement,
                           int& numNodes)
{
    int noPtsX = 3;
    int noPtsY = 4;

    ValueVector<DofType> dofVector;
    dofVector.Add(displ);

    std::vector<double> weights = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    weights[1] = (1. + 1. / sqrt(2)) / 2.;
    weights[2] = (1. + 1. / sqrt(2)) / 2.;

    std::vector<double> knotsX = {0, 0, 0, 0.5, 1, 1, 1};
    std::vector<double> knotsY = {0, 0, 0, 1, 1, 1};
    std::array<std::vector<double>, 2> knots = {{knotsX, knotsY}};

    std::array<int, 2> degree = {{2, 2}};

    std::vector<Eigen::VectorXd> controlPoints;
    controlPoints.push_back(Eigen::Vector2d({-1, 0}));
    controlPoints.push_back(Eigen::Vector2d({-1, 0.4142135623730951}));
    controlPoints.push_back(Eigen::Vector2d({-0.4142135623730951, 1}));
    controlPoints.push_back(Eigen::Vector2d({0, 1}));

    controlPoints.push_back(Eigen::Vector2d({-2.5, 0}));
    controlPoints.push_back(Eigen::Vector2d({-2.5, 0.75}));
    controlPoints.push_back(Eigen::Vector2d({-0.75, 2.5}));
    controlPoints.push_back(Eigen::Vector2d({0, 2.5}));

    controlPoints.push_back(Eigen::Vector2d({-4, 0}));
    controlPoints.push_back(Eigen::Vector2d({-4, 4}));
    controlPoints.push_back(Eigen::Vector2d({-4, 4}));
    controlPoints.push_back(Eigen::Vector2d({0, 4}));

    if (refinement == 1)
    {
        std::array<std::vector<double>, 2> knotsInserted;
        knotsInserted[0].push_back(0.25);
        knotsInserted[0].push_back(0.75);

        knotsInserted[1].push_back(0.5);

        mesh.MeshWithRefinement(coord, dofVector, knots, controlPoints, weights, degree, knotsInserted);

        numNodes = 24;
    }
    else
    {
        std::vector<NodeSimple*> controlPointsCoordinates;

        for (Eigen::VectorXd& point : controlPoints)
            controlPointsCoordinates.push_back(&mesh.mNodes.Add(point));

        std::vector<NodeSimple*> controlPointsDispl;

        for (int i = 0; i < noPtsY; i++)
            for (int j = 0; j < noPtsX; j++)
                controlPointsDispl.push_back(&mesh.mNodes.Add(Eigen::Vector2d({0., 0.})));

        Nurbs<2> nurbs_coord(knots, controlPointsCoordinates, weights, degree);
        Nurbs<2> nurbs_displ(knots, controlPointsDispl, weights, degree);

        mesh.mDofInterpolations.Insert(coord, nurbs_coord);
        mesh.mDofInterpolations.Insert(displ, nurbs_displ);

        for (int i = 2; i <= noPtsY - 1; i++)
            for (int j = 2; j <= noPtsX - 1; j++)
            {
                std::array<int, 2> knotIDsCell = {{i, j}};
                ElementIga<2> igaCoordinates(knotIDsCell, mesh.mDofInterpolations.At(coord));
                ElementIga<2> igaDispl(knotIDsCell, mesh.mDofInterpolations.At(displ));

                ElementCollectionIga<2>& elementcollection = mesh.mElements.Add(igaCoordinates);
                elementcollection.AddDofElement(displ, igaDispl);
            }

        numNodes = 12;
    }
}

BOOST_AUTO_TEST_CASE(IGA_PlateWithHoleNeumann)
{
    MeshIga<2> mesh;
    NuTo::DofType coord("coordinates", 2);
    NuTo::DofType displ("displacements", 2);
    int refinement = 1;
    int numNodes = 0;
    PlateWithHoleGeometry(mesh, coord, displ, refinement, numNodes);

    // %%%%%%%%%%%%%%%%%%%% Constraints %%%%%%%%%%%%%%%%%%%%%% //
    NuTo::Constraint::Constraints constraints;

    NuTo::Group<NuTo::NodeSimple> nodesConstrainedX = mesh.NodesAtAxis(NuTo::eDirection::X, displ, 0.);
    constraints.Add(displ, NuTo::Constraint::Component(nodesConstrainedX, {NuTo::eDirection::X}));

    NuTo::Group<NuTo::NodeSimple> nodesConstrainedY = mesh.NodesAtAxis(NuTo::eDirection::Y, displ, 0.);
    constraints.Add(displ, NuTo::Constraint::Component(nodesConstrainedY, {NuTo::eDirection::Y}));

    // %%%%%%%%%%%%%%%%%%%% All nodes aka. CPs %%%%%%%%%%%%%%%%%%%%%% //

    NuTo::Group<NuTo::NodeSimple> nodesDispl;
    for (int nodeId = 0; nodeId < numNodes; nodeId++)
        nodesDispl.Add(*mesh.mDofInterpolations.At(displ).GetControlPoint(nodeId));

    NuTo::DofInfo dofInfo = NuTo::DofNumbering::Build(nodesDispl, displ, constraints);

    // %%%%%%%%%%%%%%%%%%%%%%%% material and integrand %%%%%%%%%%%%%%%%%%%%%%%%%%%% //
    constexpr double E = 20000;
    constexpr double nu = 0.0;
    NuTo::Laws::LinearElastic<2> linearElasticLaw(E, nu, NuTo::ePlaneState::PLANE_STRESS);
    Integrands::MomentumBalance<2> momentumBalance(displ, linearElasticLaw);
    auto GradientF = Bind(momentumBalance, &Integrands::MomentumBalance<2>::Gradient);
    auto Hessian0F = Bind(momentumBalance, &Integrands::MomentumBalance<2>::Hessian0);

    boost::ptr_vector<CellInterface> cellContainer;
    int numIPs = 3;
    IntegrationTypeTensorProduct<2> integrationType(numIPs, eIntegrationMethod::GAUSS);

    Group<CellInterface> cellGroup;
    int cellId = 0;
    for (auto& element : mesh.mElements)
    {
        cellContainer.push_back(new Cell(element, integrationType, cellId++));
        cellGroup.Add(cellContainer.back());
    }

    // %%%%%%%%%%%%%%%%%%%%%%%% Neumann BC %%%%%%%%%%%%%%%%%%%%%%%%%%%% //
    ValueVector<DofType> dofVector;
    dofVector.Add(displ);
    eDirection dirNotIncluded = eDirection::Y;
    mesh.CreateBoundary(coord, dofVector, dirNotIncluded, 1);

    IntegrationTypeTensorProduct<1> integrationTypeBcLeft(3, eIntegrationMethod::GAUSS);
    Eigen::Vector2d pressureBCLeft(-1000, 0);
    Integrands::NeumannBc<2> neumannBcLeft(displ, pressureBCLeft);
    auto NeumannLoadLeft = Bind(neumannBcLeft, &Integrands::NeumannBc<2>::ExternalLoad);

    NuTo::Group<CellInterface> neumannCellsLeft;

    for (int i = 0; i < 2; i++)
    {
        ElementCollectionIga<1>& element = mesh.mElementsBoundary[i];
        cellContainer.push_back(new Cell(element, integrationTypeBcLeft, cellId++));
        neumannCellsLeft.Add(cellContainer.back());
    }

    IntegrationTypeTensorProduct<1> integrationTypeBcRight(3, eIntegrationMethod::GAUSS);
    Eigen::Vector2d pressureBCRight(0, 1000);
    Integrands::NeumannBc<2> neumannBcRight(displ, pressureBCRight);
    auto NeumannLoadRight = Bind(neumannBcRight, &Integrands::NeumannBc<2>::ExternalLoad);

    NuTo::Group<CellInterface> neumannCellsRight;

    for (int i = 2; i < 4; i++)
    {
        ElementCollectionIga<1>& element = mesh.mElementsBoundary[i];
        cellContainer.push_back(new Cell(element, integrationTypeBcRight, cellId++));
        neumannCellsRight.Add(cellContainer.back());
    }

    // %%%%%%%%%%%%%%%%%%%%%%%% assemble and solve %%%%%%%%%%%%%%%%%%%%%%%%%%%% //
    SimpleAssembler assembler(dofInfo);

    auto gradient = assembler.BuildVector(neumannCellsLeft, {displ}, NeumannLoadLeft);
    gradient += assembler.BuildVector(neumannCellsRight, {displ}, NeumannLoadRight);
    gradient += assembler.BuildVector(cellGroup, {displ}, GradientF);

    auto hessian = assembler.BuildMatrix(cellGroup, {displ}, Hessian0F);

    int numIndependentDofs = dofInfo.numIndependentDofs[displ];

    Eigen::MatrixXd CMat = constraints.BuildConstraintMatrix(displ, numIndependentDofs);

    BOOST_TEST_MESSAGE("CMat \n" << CMat);

    GlobalDofVector u = Solve(hessian, gradient, constraints, displ, numIndependentDofs, 0.0);

    // %%%%%%%%%%%%%%%%% merge dof values %%%%%%%%% //
    for (auto& node : nodesDispl)
    {
        node.SetValue(0, u(displ, node.GetDofNumber(0)));
        node.SetValue(1, u(displ, node.GetDofNumber(1)));
    }

    //    Visualize::Spacing s = Visualize::Spacing::GAUSS;
    Visualize::Visualizer visualize(cellGroup, Visualize::VoronoiHandler(Visualize::VoronoiGeometryQuad(numIPs)));
    visualize.DofValues(displ);

    auto stress = [linearElasticLaw, displ](const CellIpData& cellIpData) {
        return linearElasticLaw.Stress(cellIpData.Apply(displ, Nabla::Strain()));
    };
    visualize.CellData(stress, "Stress");

    visualize.CellData([](const CellIpData& cipd) { return EigenCompanion::ToEigen(cipd.Ids().cellId); }, "CellId");
    visualize.WriteVtuFile("plateWithHoleNeumann.vtu");
}

// BOOST_AUTO_TEST_CASE(IGA_PlateWithHoleNeumannRefinement0)
//{
//    MeshIga<2> mesh;
//    NuTo::DofType coord("coordinates", 2);
//    NuTo::DofType displ("displacements", 2);
//    int refinement = 0;
//    int numNodes = 12;
//    PlateWithHoleGeometry(mesh, coord, displ, refinement);

//    // %%%%%%%%%%%%%%%%%%%% Constraints %%%%%%%%%%%%%%%%%%%%%% //

//    NuTo::Constraint::Constraints constraints;

//    // Y symmetric boundary conditions
//    NuTo::Group<NuTo::NodeSimple> nodesConstrainedY;
//    std::vector<int> nodeIds = {0, 4, 8};

//    for (int id : nodeIds)
//        nodesConstrainedY.Add(*mesh.mDofInterpolations.At(displ).GetControlPoint(id));

//    constraints.Add(displ, NuTo::Constraint::Component(nodesConstrainedY, {NuTo::eDirection::Y}));

//    // X symmetric boundary conditions
//    NuTo::Group<NuTo::NodeSimple> nodesConstrainedX;
//    nodeIds = {3, 7, 11};

//    for (int id : nodeIds)
//        nodesConstrainedX.Add(*mesh.mDofInterpolations.At(displ).GetControlPoint(id));

//    constraints.Add(displ, NuTo::Constraint::Component(nodesConstrainedX, {NuTo::eDirection::X}));

//    // %%%%%%%%%%%%%%%%%%%% All nodes aka. CPs %%%%%%%%%%%%%%%%%%%%%% //

//    NuTo::Group<NuTo::NodeSimple> nodesDispl;
//    for (int nodeId = 0; nodeId < numNodes; nodeId++)
//        nodesDispl.Add(*mesh.mDofInterpolations.At(displ).GetControlPoint(nodeId));

//    NuTo::DofInfo dofInfo = NuTo::DofNumbering::Build(nodesDispl, displ, constraints);

//    // %%%%%%%%%%%%%%%%%%%%%%%% material and integrand %%%%%%%%%%%%%%%%%%%%%%%%%%%% //
//    constexpr double E = 20000;
//    constexpr double nu = 0.0;
//    NuTo::Laws::LinearElastic<2> linearElasticLaw(E, nu, NuTo::ePlaneState::PLANE_STRESS);
//    Integrands::MomentumBalance<2> momentumBalance(displ, linearElasticLaw);
//    auto GradientF = Bind(momentumBalance, &Integrands::MomentumBalance<2>::Gradient);
//    auto Hessian0F = Bind(momentumBalance, &Integrands::MomentumBalance<2>::Hessian0);

//    boost::ptr_vector<CellInterface> cellContainer;
//    int numIPs = 3;
//    IntegrationTypeTensorProduct<2> integrationType(numIPs, eIntegrationMethod::GAUSS);

//    Group<CellInterface> cellGroup;
//    int cellId = 0;
//    for (auto& element : mesh.mElements)
//    {
//        cellContainer.push_back(new Cell(element, integrationType, cellId++));
//        cellGroup.Add(cellContainer.back());
//    }

//    // %%%%%%%%%%%%%%%%%%%% Neumann BC %%%%%%%%%%%%%%%%%%%%%% //
//    std::vector<NodeSimple*> controlPointsCoordinatesNBC;
//    std::vector<NodeSimple*> controlPointsDisplNBC;

//    //    nodeIds = {8, 9, 10, 11};
//    nodeIds = {0, 4, 8};

//    for (int id : nodeIds)
//    {
//        controlPointsCoordinatesNBC.push_back(mesh.mDofInterpolations.At(coord).GetControlPoint(id));
//        controlPointsDisplNBC.push_back(mesh.mDofInterpolations.At(displ).GetControlPoint(id));
//    }

//    //    std::vector<double> knots1 = {0, 0, 0, 0.5, 1, 1, 1};
//    //    std::vector<double> weightsNBC = {1, 1, 1, 1};
//    std::vector<double> knots1 = {0, 0, 0, 1, 1, 1};
//    std::vector<double> weightsNBC = {1, 1, 1};
//    std::array<int, 1> degreeNBC = {{2}};

//    Nurbs<1> nurbsNBC_coord({{knots1}}, controlPointsCoordinatesNBC, weightsNBC, degreeNBC);
//    Nurbs<1> nurbsNBC_disp({{knots1}}, controlPointsDisplNBC, weightsNBC, degreeNBC);

//    std::vector<ElementCollectionIga<1>> vectorElementsNBC;
//    //    int knotIdMax = 3;
//    int knotIdMax = 3;

//    for (int i = 2; i <= knotIdMax; i++)
//    {
//        std::array<int, 1> knotIDsCell = {{i}};
//        ElementIga<1> igaCoordinates(knotIDsCell, nurbsNBC_coord);
//        ElementIga<1> igaDispl(knotIDsCell, nurbsNBC_disp);

//        vectorElementsNBC.push_back(igaCoordinates);
//        vectorElementsNBC.back().AddDofElement(displ, igaDispl);
//    }

//    IntegrationTypeTensorProduct<1> integrationTypeBc(3, eIntegrationMethod::GAUSS);
//    Eigen::Vector2d pressureBC(10, 0);
//    Integrands::NeumannBc<2> neumannBc(displ, pressureBC);
//    auto NeumannLoad = Bind(neumannBc, &Integrands::NeumannBc<2>::ExternalLoad);

//    NuTo::Group<CellInterface> neumannCells;

//    cellContainer.push_back(new Cell(vectorElementsNBC[0], integrationTypeBc, cellId++));
//    neumannCells.Add(cellContainer.back());
//    //    cellContainer.push_back(new Cell(vectorElementsNBC[1], integrationTypeBc, cellId++));
//    //    neumannCells.Add(cellContainer.back());

//    // %%%%%%%%%%%%%%%%%%%%%%%% assemble and solve %%%%%%%%%%%%%%%%%%%%%%%%%%%% //
//    SimpleAssembler assembler(dofInfo);

//    auto gradient = assembler.BuildVector(neumannCells, {displ}, NeumannLoad);
//    gradient += assembler.BuildVector(cellGroup, {displ}, GradientF);

//    auto hessian = assembler.BuildMatrix(cellGroup, {displ}, Hessian0F);

//    int numIndependentDofs = dofInfo.numIndependentDofs[displ];

//    Eigen::MatrixXd CMat = constraints.BuildConstraintMatrix(displ, numIndependentDofs);

//    BOOST_TEST_MESSAGE("CMat \n" << CMat);

//    GlobalDofVector u = Solve(hessian, gradient, constraints, displ, numIndependentDofs, 0.0);

//    // %%%%%%%%%%%%%%%%% merge dof values %%%%%%%%% //
//    for (auto& node : nodesDispl)
//    {
//        node.SetValue(0, u(displ, node.GetDofNumber(0)));
//        node.SetValue(1, u(displ, node.GetDofNumber(1)));
//    }

//    //    Visualize::Spacing s = Visualize::Spacing::GAUSS;
//    Visualize::Visualizer visualize(cellGroup, Visualize::VoronoiHandler(Visualize::VoronoiGeometryQuad(numIPs)));
//    visualize.DofValues(displ);

//    auto stress = [linearElasticLaw, displ](const CellIpData& cellIpData) {
//        return linearElasticLaw.Stress(cellIpData.Apply(displ, Nabla::Strain()));
//    };
//    visualize.CellData(stress, "Stress");

//    visualize.CellData([](const CellIpData& cipd) { return EigenCompanion::ToEigen(cipd.Ids().cellId); }, "CellId");
//    visualize.WriteVtuFile("plateWithHoleNeumannRef0.vtu");
//}
