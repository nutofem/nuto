#include "BoostUnitTest.h"
#include <boost/filesystem.hpp>
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

#include "mechanics/cell/Cell.h"
#include "mechanics/cell/SimpleAssembler.h"

#include "visualize/AverageHandler.h"
#include "visualize/AverageGeometries.h"
#include "visualize/VoronoiHandler.h"
#include "visualize/VoronoiGeometries.h"
#include "visualize/Visualizer.h"

#include "PlateWithHoleAnalytic.h"

// void ApplyBCs(NuTo::Structure& s)
//{
//    constexpr double lx = 4;
//    constexpr double ly = 4;
//
//    auto& groupLeft = s.GroupGetNodesAtCoordinate(NuTo::eDirection::X, 0.);
//    auto& groupRight = s.GroupGetNodesAtCoordinate(NuTo::eDirection::X, lx);
//    auto& groupLower = s.GroupGetNodesAtCoordinate(NuTo::eDirection::Y, 0);
//    auto& groupUpper = s.GroupGetNodesAtCoordinate(NuTo::eDirection::Y, ly);
//
//    s.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Constraint::Component(groupLeft,
//    {NuTo::eDirection::X}));
//    s.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS,
//                        NuTo::Constraint::Component(groupLower, {NuTo::eDirection::Y}));
//
//    int groupElementBCUpper = s.GroupCreate(NuTo::eGroupId::Elements);
//    int groupElementBCRight = s.GroupCreate(NuTo::eGroupId::Elements);
//    s.GroupAddElementsFromNodes(groupElementBCRight, s.GroupGetId(&groupRight), false);
//    s.GroupAddElementsFromNodes(groupElementBCUpper, s.GroupGetId(&groupUpper), false);
//
//    s.LoadSurfacePressureFunctionCreate2D(groupElementBCRight, s.GroupGetId(&groupRight),
//                                          NuTo::Test::PlateWithHoleAnalytical::PressureRight);
//    s.LoadSurfacePressureFunctionCreate2D(groupElementBCUpper, s.GroupGetId(&groupUpper),
//                                          NuTo::Test::PlateWithHoleAnalytical::PressureTop);
//}
//
// void CheckSolution(NuTo::Structure& s, double tolerance)
//{
//    s.SetShowTime(false);
//    for (int elementId : s.GroupGetMemberIds(s.GroupGetElementsTotal()))
//    {
//        auto ipCoords = s.ElementGetIntegrationPointCoordinates(elementId);
//        auto ipStress = s.ElementGetEngineeringStress(elementId);
//
//        for (int iIP = 0; iIP < ipCoords.cols(); ++iIP)
//        {
//            auto numericStressNuTo = ipStress.col(iIP);
//            Eigen::Vector3d numericStress(numericStressNuTo[0], numericStressNuTo[1], numericStressNuTo[5]);
//            auto analyticStress = NuTo::Test::PlateWithHoleAnalytical::AnalyticStress(ipCoords.col(iIP));
//            auto error = (analyticStress - numericStress).norm() / analyticStress.norm();
//            BOOST_CHECK_SMALL(error, tolerance);
//        }
//    }
//}
//

using namespace NuTo;

Constraint::Constraints FixBottomAndLeft(MeshFem* rMesh, DofType disp)
{
    Constraint::Constraints constraints;

    auto nodesLeft = rMesh->NodesAtAxis(eDirection::X, disp, 0.);
    auto nodesLower = rMesh->NodesAtAxis(eDirection::Y, disp, 0.);

    constraints.Add(disp, Constraint::Component(nodesLeft, {eDirection::X}));
    constraints.Add(disp, Constraint::Component(nodesLower, {eDirection::Y}));

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
    AddDofInterpolation(&meshGmsh.GetMeshFEM(), disp);

    auto constraints = FixBottomAndLeft(&mesh, disp);

    constexpr double lx = 4;
    constexpr double ly = 4;
    auto nodesRight = mesh.NodesAtAxis(eDirection::X, disp, lx);
    auto nodesUpper = mesh.NodesAtAxis(eDirection::Y, disp, ly);


    // auto meshInfo = s.ImportFromGmsh(meshFile);
    //
    // int interpolationType = meshInfo[0].second;
    // s.InterpolationTypeAdd(interpolationType, NuTo::Node::eDof::DISPLACEMENTS,
    //                       NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    //
    // s.SetVerboseLevel(10);
    // s.ElementTotalConvertToInterpolationType();
    //
    // double thickness = 1.;
    // auto section = NuTo::SectionPlane::Create(thickness, false);
    // s.ElementTotalSetSection(section);
    //
    // using namespace NuTo::Constitutive;
    //
    // int constitutiveLaw = s.ConstitutiveLawCreate(eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    // s.ConstitutiveLawSetParameterDouble(constitutiveLaw, eConstitutiveParameter::YOUNGS_MODULUS, 1.e5);
    // s.ConstitutiveLawSetParameterDouble(constitutiveLaw, eConstitutiveParameter::POISSONS_RATIO, 0.3);
    // s.ElementTotalSetConstitutiveLaw(constitutiveLaw);
    //
    // ApplyBCs(s);
    //
    // s.SolveGlobalSystemStaticElastic();
    //
    // int visualizationGroup = s.GroupGetElementsTotal();
    // s.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::DISPLACEMENTS);
    // s.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);
    // s.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRESS);
    //
    // std::string resultDir = "./PlateWithHoleResults";
    // boost::filesystem::create_directory(resultDir);
    // s.ExportVtkDataFileElements(resultDir + "/PlateWithHole.vtu", true);
    //
    // CheckSolution(s, 0.05);
}
