#include "BoostUnitTest.h"
#include <string>
#include <iostream>
#include <boost/filesystem.hpp>

#include "mechanics/structures/unstructured/Structure.h"

#include "mechanics/groups/Group.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/MechanicsEnums.h"
#include "visualize/VisualizeEnum.h"
#include "mechanics/sections/SectionPlane.h"
#include "mechanics/constraints/ConstraintCompanion.h"

#include "PlateWithHoleAnalytic.h"

void ApplyBCs(NuTo::Structure& s)
{
    constexpr double lx        = 4;
    constexpr double ly        = 4;

    auto& groupLeft  = s.GroupGetNodesAtCoordinate(NuTo::eDirection::X, 0.);
    auto& groupRight = s.GroupGetNodesAtCoordinate(NuTo::eDirection::X, lx);
    auto& groupLower = s.GroupGetNodesAtCoordinate(NuTo::eDirection::Y, 0);
    auto& groupUpper = s.GroupGetNodesAtCoordinate(NuTo::eDirection::Y, ly);

    s.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS,
            NuTo::Constraint::Component(groupLeft, {NuTo::eDirection::X}));
    s.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS,
            NuTo::Constraint::Component(groupLower, {NuTo::eDirection::Y}));

    int groupElementBCUpper = s.GroupCreate(NuTo::eGroupId::Elements);
    int groupElementBCRight = s.GroupCreate(NuTo::eGroupId::Elements);
    s.GroupAddElementsFromNodes(groupElementBCRight, s.GroupGetId(&groupRight), false);
    s.GroupAddElementsFromNodes(groupElementBCUpper, s.GroupGetId(&groupUpper), false);
    
    s.LoadSurfacePressureFunctionCreate2D(groupElementBCRight, s.GroupGetId(&groupRight), NuTo::Test::PlateWithHoleAnalytical::PressureRight);
    s.LoadSurfacePressureFunctionCreate2D(groupElementBCUpper, s.GroupGetId(&groupUpper), NuTo::Test::PlateWithHoleAnalytical::PressureTop);
}

void CheckSolution(NuTo::Structure& s, double tolerance)
{
    s.SetShowTime(false);
    for (int elementId : s.GroupGetMemberIds(s.GroupGetElementsTotal()))
    {
        auto ipCoords = s.ElementGetIntegrationPointCoordinates(elementId);
        auto ipStress = s.ElementGetEngineeringStress(elementId);

        for (int iIP = 0; iIP < ipCoords.cols(); ++iIP)
        {
            auto numericStressNuTo = ipStress.col(iIP);
            Eigen::Vector3d numericStress(numericStressNuTo[0], numericStressNuTo[1], numericStressNuTo[5]);
            auto analyticStress  = NuTo::Test::PlateWithHoleAnalytical::AnalyticStress(ipCoords.col(iIP));
            auto error = (analyticStress - numericStress).norm() / analyticStress.norm();
            BOOST_CHECK_SMALL(error, tolerance);
        }
    }
}

BOOST_AUTO_TEST_CASE(PlateWithHole)
{
    auto binary = boost::unit_test::framework::master_test_suite().argv[0];
    boost::filesystem::path binaryPath = std::string(binary);
    binaryPath.remove_filename();

    std::string meshFile = binaryPath.string() + "/meshes/PlateWithHole.msh";

    NuTo::Structure s(2);
    auto meshInfo = s.ImportFromGmsh(meshFile);

    int interpolationType = meshInfo[0].second;
    s.InterpolationTypeAdd(interpolationType, NuTo::Node::eDof::DISPLACEMENTS,
                           NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);

    s.SetVerboseLevel(10);
    s.ElementTotalConvertToInterpolationType();

    double thickness = 1.;
    auto section = NuTo::SectionPlane::Create(thickness, false);
    s.ElementTotalSetSection(section);

    using namespace NuTo::Constitutive;

    int constitutiveLaw = s.ConstitutiveLawCreate(eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    s.ConstitutiveLawSetParameterDouble(constitutiveLaw, eConstitutiveParameter::YOUNGS_MODULUS, 1.e5);
    s.ConstitutiveLawSetParameterDouble(constitutiveLaw, eConstitutiveParameter::POISSONS_RATIO, 0.3);
    s.ElementTotalSetConstitutiveLaw(constitutiveLaw);

    ApplyBCs(s);

    s.SolveGlobalSystemStaticElastic();

    int visualizationGroup = s.GroupGetElementsTotal();
    s.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::DISPLACEMENTS);
    s.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);
    s.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRESS);

    std::string resultDir = "./PlateWithHoleResults";
    boost::filesystem::create_directory(resultDir);
    s.ExportVtkDataFileElements(resultDir + "/PlateWithHole.vtu", true);

    CheckSolution(s, 0.05);
}
