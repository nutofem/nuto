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
#include "nuto/math/EigenCompanion.h"

#include "nuto/mechanics/mesh/MeshIga.h"
using namespace NuTo;

/*  |>*----*----* --> *
      |    |    |
      |    |    |
    |>*----*----* --> *
      |    |    |
      |    |    |
    |>*----*----* --> *
      /\
      --
*/

void QuadPatchTestNurbs(MeshIga<2>& mesh, const NuTo::DofType& coord, const NuTo::DofType& displ)
{
    std::vector<NodeSimple*> controlPointsCoordinates;

    int num = 3;
    for (int i = 0; i < num; i++)
        for (int j = 0; j < num; j++)
        {
            auto& node = mesh.mNodes.Add(Eigen::Vector2d({j, i}));
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

BOOST_AUTO_TEST_CASE(PatchTestIgaMesh)
{
    MeshIga<2> mesh;
    NuTo::DofType coord("coordinates", 2);
    NuTo::DofType displ("displacements", 2);
    QuadPatchTestNurbs(mesh, coord, displ);

    // %%%%%%%%%%%%%%%%%%%% Constraints %%%%%%%%%%%%%%%%%%%%%% //
    NuTo::Constraint::Constraints constraints;

    Eigen::Vector2i nodeIds;
    nodeIds << 0, 0;
    NuTo::Group<NuTo::NodeSimple> nodesConstrained;
    nodesConstrained.Add(*mesh.mDofInterpolations.At(displ).GetControlPoint(nodeIds));
    nodeIds << 0, 1;
    nodesConstrained.Add(*mesh.mDofInterpolations.At(displ).GetControlPoint(nodeIds));
    nodeIds << 0, 2;
    nodesConstrained.Add(*mesh.mDofInterpolations.At(displ).GetControlPoint(nodeIds));
    constraints.Add(displ, NuTo::Constraint::Component(nodesConstrained, {NuTo::eDirection::X}));

    NuTo::Group<NuTo::NodeSimple> nodesConstrained1;
    nodeIds << 0, 0;
    nodesConstrained1.Add(*mesh.mDofInterpolations.At(displ).GetControlPoint(nodeIds));
    constraints.Add(displ, NuTo::Constraint::Component(nodesConstrained1, {NuTo::eDirection::Y}));

    NuTo::Group<NuTo::NodeSimple> nodesDisplConstrained;
    nodeIds << 2, 0;
    nodesDisplConstrained.Add(*mesh.mDofInterpolations.At(displ).GetControlPoint(nodeIds));
    nodeIds << 2, 1;
    nodesDisplConstrained.Add(*mesh.mDofInterpolations.At(displ).GetControlPoint(nodeIds));
    nodeIds << 2, 2;
    nodesDisplConstrained.Add(*mesh.mDofInterpolations.At(displ).GetControlPoint(nodeIds));
    constraints.Add(displ, NuTo::Constraint::Component(nodesDisplConstrained, {NuTo::eDirection::X}, 1));

    NuTo::Group<NuTo::NodeSimple> nodesDispl;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        {
            nodeIds << i, j;
            nodesDispl.Add(*mesh.mDofInterpolations.At(displ).GetControlPoint(nodeIds));
        }


    NuTo::DofInfo dofInfo = NuTo::DofNumbering::Build(nodesDispl, displ, constraints);
    const int numDofs = dofInfo.numIndependentDofs[displ] + dofInfo.numDependentDofs[displ];
    const int numDepDofs = dofInfo.numDependentDofs[displ];
    Eigen::MatrixXd CMat = constraints.BuildConstraintMatrix(displ, numDofs - numDepDofs);

    BOOST_TEST_MESSAGE("CMat \n" << CMat);

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
    cellContainer.push_back(new Cell(mesh.mElements[0], integrationType, cellId++));
    cellGroup.Add(cellContainer.back());

    // %%%%%%%%%%%%%%%%%%%%%%%% assemble and solve %%%%%%%%%%%%%%%%%%%%%%%%%%%% //
    SimpleAssembler assembler(dofInfo);

    auto gradient = assembler.BuildVector(cellGroup, {displ}, GradientF);
    auto hessian = assembler.BuildMatrix(cellGroup, {displ}, Hessian0F);

    int numIndependentDofs = dofInfo.numIndependentDofs[displ];
    GlobalDofVector u = Solve(hessian, gradient, constraints, displ, numIndependentDofs, 0.0);

    // %%%%%%%%%%%%%%%%% merge dof values %%%%%%%%% //
    for (auto& node : nodesDispl)
    {
        node.SetValue(0, u(displ, node.GetDofNumber(0)));
        node.SetValue(1, u(displ, node.GetDofNumber(1)));
    }

    Visualize::Visualizer visualize(cellGroup, Visualize::VoronoiHandler(Visualize::VoronoiGeometryQuad(numIPs)));
    visualize.DofValues(displ);

    auto stress = [linearElasticLaw, displ](const CellIpData& cellIpData) {
        return linearElasticLaw.Stress(cellIpData.Apply(displ, Nabla::Strain()));
    };
    visualize.CellData(stress, "Stress");

    visualize.CellData([](const CellIpData&) { return EigenCompanion::ToEigen(7.0); }, "Seven");
    visualize.CellData([](const CellIpData& cipd) { return EigenCompanion::ToEigen(cipd.Ids().cellId); }, "CellId");
    visualize.WriteVtuFile("outputVoronoi.vtu");
}

// BOOST_AUTO_TEST_CASE(IGA_PlateWithHoleNeumann)
//{
//    NuTo::Structure* s = new NuTo::Structure(2);

//    /*********/
//    // Mesh  //
//    /*********/

//    int refine = 5;
//    int noPtsX = 3;
//    int noPtsY = 4;

//    Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic> controlPts(noPtsX, noPtsY);

//    controlPts(0, 0) = Eigen::Vector2d(-1, 0);
//    controlPts(0, 1) = Eigen::Vector2d(-1, 0.4142135623730951);
//    controlPts(0, 2) = Eigen::Vector2d(-0.4142135623730951, 1);
//    controlPts(0, 3) = Eigen::Vector2d(0, 1);

//    controlPts(1, 0) = Eigen::Vector2d(-2.5, 0);
//    controlPts(1, 1) = Eigen::Vector2d(-2.5, 0.75);
//    controlPts(1, 2) = Eigen::Vector2d(-0.75, 2.5);
//    controlPts(1, 3) = Eigen::Vector2d(0, 2.5);

//    controlPts(2, 0) = Eigen::Vector2d(-4, 0);
//    controlPts(2, 1) = Eigen::Vector2d(-4, 4);
//    controlPts(2, 2) = Eigen::Vector2d(-4, 4);
//    controlPts(2, 3) = Eigen::Vector2d(0, 4);

//    Eigen::MatrixXd weights(noPtsX, noPtsY);
//    weights.setOnes(noPtsX, noPtsY);
//    weights(0, 1) = (1. + 1. / sqrt(2)) / 2.;
//    weights(0, 2) = (1. + 1. / sqrt(2)) / 2.;


//    Eigen::VectorXd knotsX(7);
//    knotsX << 0, 0, 0, 0.5, 1, 1, 1;
//    Eigen::VectorXd knotsY(6);
//    knotsY << 0, 0, 0, 1, 1, 1;

//    Eigen::Vector2i degree(2, 2);

//    NuTo::BSplineSurface surface(degree, knotsX, knotsY, controlPts, weights);

//    for (int i = 0; i < refine; i++)
//    {
//        surface.DuplicateKnots(0);
//        surface.DuplicateKnots(1);
//    }

//    /**************/
//    // Structure  //
//    /**************/

//    std::set<NuTo::Node::eDof> setOfDOFS;
//    setOfDOFS.insert(NuTo::Node::eDof::COORDINATES);
//    setOfDOFS.insert(NuTo::Node::eDof::DISPLACEMENTS);

//    int groupNodes = s->GroupCreate("Nodes");
//    int groupElements = s->GroupCreate("Elements");

//    surface.buildIGAStructure(*s, setOfDOFS, groupElements, groupNodes);

//    /** create section **/
//    double thickness = 1.;
//    auto mySection = NuTo::SectionPlane::Create(thickness, false);

//    /** create constitutive law **/
//    double YoungsModulus = 1.e5;
//    double PoissonRatio = 0.3;
//    int myMatLin = s->ConstitutiveLawCreate("LINEAR_ELASTIC_ENGINEERING_STRESS");
//    s->ConstitutiveLawSetParameterDouble(myMatLin, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS,
//                                         YoungsModulus);
//    s->ConstitutiveLawSetParameterDouble(myMatLin, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO,
//                                         PoissonRatio);

//    s->ElementTotalSetConstitutiveLaw(myMatLin);
//    s->ElementTotalSetSection(mySection);

//    /**********************/
//    // Boundary condition //
//    /**********************/

//    int groupElementsLeft = s->GroupCreate("ELEMENTS");
//    int groupElementsUpper = s->GroupCreate("ELEMENTS");

//    int start = (surface.GetNumIGAElements(1) - 1) * surface.GetNumIGAElements(0);
//    for (int i = start; i < start + surface.GetNumIGAElements(0) / 2; i++)
//    {
//        s->GroupAddElement(groupElementsLeft, i);
//    }

//    start += surface.GetNumIGAElements(0) / 2;
//    for (int i = start; i < surface.GetNumIGAElements(); i++)
//    {
//        s->GroupAddElement(groupElementsUpper, i);
//    }

//    auto& groupLeft = s->GroupGetNodesAtCoordinate(NuTo::eDirection::X, -4);
//    auto& groupRight = s->GroupGetNodesAtCoordinate(NuTo::eDirection::X, 0);
//    auto& groupBottom = s->GroupGetNodesAtCoordinate(NuTo::eDirection::Y, 0);
//    auto& groupTop = s->GroupGetNodesAtCoordinate(NuTo::eDirection::Y, 4);

//    s->Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS,
//                         NuTo::Constraint::Component(groupBottom, {NuTo::eDirection::Y}));
//    s->Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS,
//                         NuTo::Constraint::Component(groupRight, {NuTo::eDirection::X}));

//    s->LoadSurfacePressureFunctionCreate2D(groupElementsLeft, s->GroupGetId(&groupLeft),
//                                           NuTo::Test::PlateWithHoleAnalytical::PressureLeft);
//    s->LoadSurfacePressureFunctionCreate2D(groupElementsUpper, s->GroupGetId(&groupTop),
//                                           NuTo::Test::PlateWithHoleAnalytical::PressureTop);

//    BOOST_CHECK_NO_THROW(SolveAndVisualize(s, "Hole5"));
//}
