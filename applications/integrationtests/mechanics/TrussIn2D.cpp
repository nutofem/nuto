/*
 *  Test
 *
 *  Created on: 11 August 2015
 *      Author: phuschke
 *
 *
 *  Compares the results of three identical trusses
 *  in 2D at different inclinations (0°, 45°, and 90°).
 *
 *
 */
#include <boost/filesystem.hpp>
#include "mechanics/MechanicsEnums.h"
#include "mechanics/constraints/ConstraintCompanion.h"
#include "mechanics/sections/SectionTruss.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"


class ParametersMaterial
{
public:
    static constexpr double mYoungsModulus = 10;
    static constexpr double mPoissonsRatio = 0.0;
};

class ParametersTimeIntegration
{
public:
    static constexpr double mTimeStep = 1.;
    static constexpr double mSimulationTime = 1.0;
    static constexpr double mLoad = 1.0;
};

class ParametersGeometry2D
{
public:
    static const int mDimension;
    static const double mCrossSection;
    static const Eigen::Vector2d mDirectionX;
    static const Eigen::Vector2d mDirectionY;
};

const int ParametersGeometry2D::mDimension = 2;
const double ParametersGeometry2D::mCrossSection = 0.1;
const Eigen::Vector2d ParametersGeometry2D::mDirectionX = Eigen::Vector2d::UnitX();
const Eigen::Vector2d ParametersGeometry2D::mDirectionY = Eigen::Vector2d::UnitY();

void Run2d(Eigen::VectorXd rNodeCoords0, Eigen::VectorXd rNodeCoords1, Eigen::VectorXd rNodeCoords2,
           Eigen::VectorXd rDirectionAligned, Eigen::VectorXd rDirectionOrthogonal)
{
    //**********************************************
    //          Structure
    //**********************************************

    NuTo::Structure myStructure(ParametersGeometry2D::mDimension);
    myStructure.SetVerboseLevel(10);

    //**********************************************
    //         Integration Scheme
    //**********************************************

    boost::filesystem::path resultPath(boost::filesystem::initial_path().string() + "/resultTrussIn2D/");
    boost::filesystem::remove_all(resultPath);
    boost::filesystem::create_directory(resultPath);

    NuTo::NewmarkDirect myIntegrationScheme(&myStructure);
    myIntegrationScheme.SetTimeStep(ParametersTimeIntegration::mTimeStep);
    myIntegrationScheme.PostProcessing().SetResultDirectory(resultPath.string(), false);

    //**********************************************
    //          Section
    //**********************************************

    auto fibreSection = NuTo::SectionTruss::Create(ParametersGeometry2D::mCrossSection);

    //**********************************************
    //          Material
    //**********************************************

    int fibreMaterial =
            myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    myStructure.ConstitutiveLawSetParameterDouble(fibreMaterial,
                                                  NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS,
                                                  ParametersMaterial::mYoungsModulus);
    myStructure.ConstitutiveLawSetParameterDouble(fibreMaterial,
                                                  NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO,
                                                  ParametersMaterial::mPoissonsRatio);

    //**********************************************
    //          Interpolation
    //**********************************************

    int fibreInterpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TRUSSXD);
    myStructure.InterpolationTypeAdd(fibreInterpolationType, NuTo::Node::eDof::COORDINATES,
                                     NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    myStructure.InterpolationTypeAdd(fibreInterpolationType, NuTo::Node::eDof::DISPLACEMENTS,
                                     NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    myStructure.InterpolationTypeSetIntegrationType(fibreInterpolationType,
                                                    NuTo::eIntegrationType::IntegrationType1D2NGauss2Ip);

    //**********************************************
    //          Geometry
    //**********************************************

    // Nodes
    int nodeId0 = myStructure.NodeCreate(rNodeCoords0);
    int nodeId1 = myStructure.NodeCreate(rNodeCoords1);
    int nodeId2 = myStructure.NodeCreate(rNodeCoords2);

    std::vector<int> nodeIndices({nodeId0, nodeId1, nodeId2});
    myStructure.ElementCreate(fibreInterpolationType, nodeIndices);

    myStructure.ElementTotalConvertToInterpolationType(1e-6, 10);
    myStructure.ElementTotalSetSection(fibreSection);
    myStructure.ElementTotalSetConstitutiveLaw(fibreMaterial);

    //**********************************************
    //          Boundary Conditions
    //**********************************************
    const auto& node0 = *myStructure.NodeGetNodePtr(nodeId0);
    const auto& node1 = *myStructure.NodeGetNodePtr(nodeId1);
    const auto& node2 = *myStructure.NodeGetNodePtr(nodeId2);
    myStructure.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS,
                                  NuTo::Constraint::Component(node0, {NuTo::eDirection::X, NuTo::eDirection::Y}));
    myStructure.Constraints().Add(
            NuTo::Node::eDof::DISPLACEMENTS,
            NuTo::Constraint::Direction(node1, rDirectionOrthogonal, NuTo::Constraint::RhsConstant(0)));
    myStructure.Constraints().Add(
            NuTo::Node::eDof::DISPLACEMENTS,
            NuTo::Constraint::Direction(node2, rDirectionOrthogonal, NuTo::Constraint::RhsConstant(0)));

    //**********************************************
    //          Loads
    //**********************************************

    int load = myStructure.LoadCreateNodeForce(nodeId2, rDirectionAligned, 1);

    Eigen::Matrix2d timeDependentLoad;
    timeDependentLoad(0, 0) = 0;
    timeDependentLoad(1, 0) = ParametersTimeIntegration::mSimulationTime;
    timeDependentLoad(0, 1) = 0;
    timeDependentLoad(1, 1) = ParametersTimeIntegration::mLoad;

    myIntegrationScheme.SetTimeDependentLoadCase(load, timeDependentLoad);

    //**********************************************
    //          Solver
    //**********************************************

    myStructure.NodeBuildGlobalDofs();
    myStructure.CalculateMaximumIndependentSets();

    myIntegrationScheme.Solve(ParametersTimeIntegration::mSimulationTime);

    Eigen::VectorXd displacements;
    boost::filesystem::remove_all(resultPath);

    myStructure.NodeGetDisplacements(nodeId2, displacements);

    if (std::abs(displacements.norm() - std::sqrt(2)) > 2e-8)
    {
        std::cout << myStructure.Constraints() << std::endl;
        std::cout << "Direction: " << rDirectionAligned.transpose() << std::endl;
        std::cout << "Displacements: " << displacements.transpose() << std::endl;
        throw NuTo::Exception("The calculated displacements do not agree with the analytical solution!");
    }
}

int main(int argc, char* argv[])
{

    // node coordinates
    Eigen::VectorXd nodeCoords0(ParametersGeometry2D::mDimension);
    Eigen::VectorXd nodeCoords1(ParametersGeometry2D::mDimension);
    Eigen::VectorXd nodeCoords2(ParametersGeometry2D::mDimension);

    // directions
    Eigen::VectorXd directionAligned(ParametersGeometry2D::mDimension);
    Eigen::VectorXd directionOrthogonal(ParametersGeometry2D::mDimension);

    try
    {
        std::cout << "////////////////////////////////////////////////" << std::endl;
        std::cout << "example 1: quadratic truss element, 0° inclined" << std::endl;
        std::cout << "////////////////////////////////////////////////" << std::endl;


        nodeCoords0[0] = 1.0;
        nodeCoords0[1] = 2.0;

        nodeCoords1[0] = 1.0 + 0.5 * std::sqrt(2);
        nodeCoords1[1] = 2.0;

        nodeCoords2[0] = 1.0 + 1.0 * std::sqrt(2);
        nodeCoords2[1] = 2.0;

        directionAligned[0] = 1.0;
        directionAligned[1] = 0.0;

        directionOrthogonal[0] = 0.0;
        directionOrthogonal[1] = 1.0;

        Run2d(nodeCoords0, nodeCoords1, nodeCoords2, directionAligned, directionOrthogonal);


        std::cout << "////////////////////////////////////////////////" << std::endl;
        std::cout << "example 2: quadratic truss element, 45° inclined" << std::endl;
        std::cout << "////////////////////////////////////////////////" << std::endl;

        nodeCoords0[0] = 1.0;
        nodeCoords0[1] = 2.0;

        nodeCoords1[0] = 1.0 + 0.5;
        nodeCoords1[1] = 2.0 + 0.5;

        nodeCoords2[0] = 1.0 + 1.0;
        nodeCoords2[1] = 2.0 + 1.0;

        directionAligned[0] = 1.0;
        directionAligned[1] = 1.0;

        directionOrthogonal[0] = -1.0;
        directionOrthogonal[1] = 1.0;

        Run2d(nodeCoords0, nodeCoords1, nodeCoords2, directionAligned, directionOrthogonal);


        std::cout << "////////////////////////////////////////////////" << std::endl;
        std::cout << "example 3: quadratic truss element, 90° inclined" << std::endl;
        std::cout << "////////////////////////////////////////////////" << std::endl;


        nodeCoords0[0] = 1.0;
        nodeCoords0[1] = 2.0;

        nodeCoords1[0] = 1.0;
        nodeCoords1[1] = 2.0 + 0.5 * std::sqrt(2);

        nodeCoords2[0] = 1.0;
        nodeCoords2[1] = 2.0 + 1.0 * std::sqrt(2);

        directionAligned[0] = 0.0;
        directionAligned[1] = 1.0;

        directionOrthogonal[0] = 1.0;
        directionOrthogonal[1] = 0.0;

        Run2d(nodeCoords0, nodeCoords1, nodeCoords2, directionAligned, directionOrthogonal);
    }
    catch (NuTo::Exception& e)
    {
        std::cout << "## Test failed ##" << std::endl;
        std::cout << e.ErrorMessage();
        return EXIT_FAILURE;
    }

    std::cout << "## Test successful ##" << std::endl;
}
