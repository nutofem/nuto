#include <map>
#include <string>
#include <iostream>
#include <cmath>
#include <boost/filesystem.hpp>

#include "mechanics/structures/unstructured/Structure.h"

#include "mechanics/nodes/NodeBase.h"
#include "mechanics/MechanicsEnums.h"
#include "visualize/VisualizeEnum.h"
#include "mechanics/sections/SectionPlane.h"
#include "mechanics/constraints/ConstraintCompanion.h"


using ExactStress     = Eigen::Vector3d;
constexpr double load = 10;


// @brief returns the analytical solution (sxx, syy, sxy) taken from
// https://en.wikiversity.org/wiki/Introduction_to_Elasticity/Plate_with_hole_in_tension
ExactStress GetAnalyticSolution(Eigen::VectorXd rCartesianCoordinate)
{
    constexpr double a = 1.0;
    double r           = rCartesianCoordinate.norm();
    double theta       = std::atan(rCartesianCoordinate[1] / rCartesianCoordinate[0]);

    double cos2t = std::cos(2 * theta);
    double cos4t = std::cos(4 * theta);
    double sin2t = std::sin(2 * theta);
    double sin4t = std::sin(4 * theta);

    double fac1 = (a * a) / (r * r);
    double fac2 = 1.5 * fac1 * fac1;
    ExactStress stress;
    stress[0] = 1. - fac1 * (1.5 * cos2t + cos4t) + fac2 * cos4t;
    stress[1] = -fac1 * (0.5 * cos2t - cos4t) - fac2 * cos4t;
    stress[2] = -fac1 * (0.5 * sin2t + sin4t) + fac2 * sin4t;

    return stress * load;
}


Eigen::Vector2d GetPressure(Eigen::VectorXd rCartesianCoordinate, Eigen::Vector2d rN)
{
    ExactStress s = GetAnalyticSolution(rCartesianCoordinate);
    Eigen::Matrix2d stress;
    stress << s[0], s[2], s[2], s[1];
    Eigen::Vector2d pressure = stress * rN;
    std::cout << "Apply pressure " << pressure.transpose() << " at coordinate " << rCartesianCoordinate.transpose()
              << std::endl;
    return pressure;
}

Eigen::Vector2d GetPressureRight(Eigen::VectorXd rCartesianCoordinate)
{
    return GetPressure(rCartesianCoordinate, Eigen::Vector2d::UnitX());
}

Eigen::Vector2d GetPressureUpper(Eigen::VectorXd rCartesianCoordinate)
{
    return GetPressure(rCartesianCoordinate, Eigen::Vector2d::UnitY());
}

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
    
    s.SetNumLoadCases(1);
    s.LoadSurfacePressureFunctionCreate2D(0, groupElementBCRight, s.GroupGetId(&groupRight), GetPressureRight);
    s.LoadSurfacePressureFunctionCreate2D(0, groupElementBCUpper, s.GroupGetId(&groupUpper), GetPressureUpper);
}

bool CheckSolution(NuTo::Structure& s)
{
    s.SetShowTime(false);
    double maxerror = 0;
    bool error = false;
    for (int elementId : s.GroupGetMemberIds(s.GroupGetElementsTotal()))
    {
        auto ipCoords = s.ElementGetIntegrationPointCoordinates(elementId);
        auto ipStress = s.ElementGetEngineeringStress(elementId);

        for (int iIP = 0; iIP < ipCoords.cols(); ++iIP)
        {
            auto numericStressNuTo = ipStress.col(iIP);
            Eigen::Vector3d numericStress(numericStressNuTo[0], numericStressNuTo[1], numericStressNuTo[5]);
            auto analyticStress  = GetAnalyticSolution(ipCoords.col(iIP));
            double error         = (numericStress - analyticStress).norm();
            double relativeError = error / analyticStress.norm();
            if (relativeError > 0.05)
            {
                std::cout << "relative error " << relativeError << '\t';
                std::cout << "N " << numericStress.transpose() << '\t';
                std::cout << "A " << analyticStress.transpose() << std::endl;
                error = true;
            }
            maxerror = std::max(maxerror, relativeError);
        }
    }
    std::cout << "maxerror " << maxerror << std::endl;
    return error;
}

//! @brief Imports a mesh file and builds the hessian and the internal gradient.
int main(int argc, char* argv[])
{
    boost::filesystem::path binaryPath = std::string(argv[0]);
    binaryPath.remove_filename();

    std::string meshFile = binaryPath.string() + "/PlateWithHole.msh";

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

    return CheckSolution(s);
}
