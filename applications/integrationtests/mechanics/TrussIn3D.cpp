/*
 *  Test
 *
 *  Created on: 11 August 2015
 *      Author: phuschke
 *
 *
 *  Compares the results of four identical trusses
 *  in 3D at different inclinations.
 *
 *
 */
#include "BoostUnitTest.h"

#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/integrationtypes/IntegrationTypeEnum.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/sections/SectionTruss.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "mechanics/constraints/ConstraintCompanion.h"

#include <boost/filesystem.hpp>
#include <iostream>
#include <fstream>
#include <string>

//**********************************************
//          Parameters
//**********************************************

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

class ParametersGeometry3D
{
public:
    static constexpr int mDimension = 3;
    static constexpr double mCrossSection = 0.1;
};



// Input:   3 node vectors (quadratic truss element),
//          1 direction vector that points in the direction of the truss (for the load)
//          2 orthogonal direction vectors (for the BC)
void Run3d(Eigen::VectorXd rNodeCoords0,
           Eigen::VectorXd rNodeCoords1,
           Eigen::VectorXd rNodeCoords2,
           Eigen::VectorXd rDirectionAligned,
           Eigen::VectorXd rDirectionOrthogonal0,
           Eigen::VectorXd rDirectionOrthogonal1,
           double expectedNorm)
{
    //**********************************************
    //          Structure
    //**********************************************

    NuTo::Structure myStructure(ParametersGeometry3D::mDimension);

    //**********************************************
    //         Integration Scheme
    //**********************************************

    boost::filesystem::path resultPath(boost::filesystem::initial_path().string() + "/resultTrussIn3D/");
    boost::filesystem::remove_all(resultPath);
    boost::filesystem::create_directory(resultPath);

    NuTo::NewmarkDirect myIntegrationScheme(&myStructure);
    myIntegrationScheme.SetTimeStep(ParametersTimeIntegration::mTimeStep);
    myIntegrationScheme.SetResultDirectory(resultPath.string(), false);

    //**********************************************
    //          Section
    //**********************************************

    auto fibreSection = NuTo::SectionTruss::Create(ParametersGeometry3D::mCrossSection);

    //**********************************************
    //          Material
    //**********************************************

    int fibreMaterial = myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    myStructure.ConstitutiveLawSetParameterDouble(fibreMaterial, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, ParametersMaterial::mYoungsModulus);
    myStructure.ConstitutiveLawSetParameterDouble(fibreMaterial, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, ParametersMaterial::mPoissonsRatio);

    //**********************************************
    //          Interpolation
    //**********************************************

    int fibreInterpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TRUSSXD);
    myStructure.InterpolationTypeAdd(fibreInterpolationType, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    myStructure.InterpolationTypeAdd(fibreInterpolationType, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    myStructure.InterpolationTypeSetIntegrationType(fibreInterpolationType, NuTo::eIntegrationType::IntegrationType1D2NGauss2Ip);

    //**********************************************
    //          Geometry
    //**********************************************

    // Nodes
    int node0 = myStructure.NodeCreate(rNodeCoords0);
    int node1 = myStructure.NodeCreate(rNodeCoords1);
    int node2 = myStructure.NodeCreate(rNodeCoords2);

    std::vector<int> nodeIndices({node0, node1, node2});
    myStructure.ElementCreate(fibreInterpolationType, nodeIndices);

    myStructure.ElementTotalConvertToInterpolationType(1e-6, 10);
    myStructure.ElementTotalSetSection(fibreSection);
    myStructure.ElementTotalSetConstitutiveLaw(fibreMaterial);

    //**********************************************
    //          Boundary Conditions
    //**********************************************
    const auto& node0Ref = *myStructure.NodeGetNodePtr(node0);
    const auto& node1Ref = *myStructure.NodeGetNodePtr(node1);
    const auto& node2Ref = *myStructure.NodeGetNodePtr(node2);
    myStructure.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS,
            NuTo::Constraint::Component(node0Ref, {NuTo::eDirection::X, NuTo::eDirection::Y, NuTo::eDirection::Z}));

    
    myStructure.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS,
            NuTo::Constraint::Direction(node1Ref, rDirectionOrthogonal0)); 
    myStructure.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS,
            NuTo::Constraint::Direction(node2Ref, rDirectionOrthogonal0)); 
    
    myStructure.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS,
            NuTo::Constraint::Direction(node1Ref, rDirectionOrthogonal1)); 
    myStructure.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS,
            NuTo::Constraint::Direction(node2Ref, rDirectionOrthogonal1)); 
   
    std::cout << myStructure.Constraints() << std::endl;

    //**********************************************
    //          Loads
    //**********************************************

    int load = myStructure.LoadCreateNodeForce(0, node2, rDirectionAligned, 1);

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

    myStructure.NodeGetDisplacements(node2, displacements);
    
    BOOST_CHECK_CLOSE(displacements.norm(), expectedNorm, 1e-8);
}

BOOST_AUTO_TEST_CASE(TrussAlignedX)
{
    std::cout << "*******************************************************************" << std::endl;
    std::cout << "*********** Start Example 1: Truss aligned with x-axis ************" << std::endl;
    std::cout << "*******************************************************************" << std::endl;

    Eigen::Vector3d directionAligned = Eigen::Vector3d::UnitX();
    Eigen::Vector3d directionOrthogonal0 = Eigen::Vector3d::UnitY();
    Eigen::Vector3d directionOrthogonal1 = Eigen::Vector3d::UnitZ();

    Run3d(Eigen::Vector3d(1, 5, 3), 
          Eigen::Vector3d(2, 5, 3), 
          Eigen::Vector3d(3, 5, 3), directionAligned, directionOrthogonal0, directionOrthogonal1, 2);

    std::cout << "*******************************************************************" << std::endl;
}

BOOST_AUTO_TEST_CASE(TrussAlignedY)
{

    std::cout << "*******************************************************************" << std::endl;
    std::cout << "*********** Start Example 2: Truss aligned with y-axis ************" << std::endl;
    std::cout << "*******************************************************************" << std::endl;
 
    Eigen::Vector3d directionAligned = Eigen::Vector3d::UnitY();
    Eigen::Vector3d directionOrthogonal0 = Eigen::Vector3d::UnitX();
    Eigen::Vector3d directionOrthogonal1 = Eigen::Vector3d::UnitZ();
 
    Run3d(Eigen::Vector3d(0,0,0), 
          Eigen::Vector3d(0,1,0), 
          Eigen::Vector3d(0,2,0), directionAligned, directionOrthogonal0, directionOrthogonal1, 2);

}

BOOST_AUTO_TEST_CASE(TrussAlignedZ)
{
    std::cout << "*******************************************************************" << std::endl;
    std::cout << "*********** Start Example 3: Truss aligned with z-axis ************" << std::endl;
    std::cout << "*******************************************************************" << std::endl;
 
    Eigen::Vector3d directionAligned = Eigen::Vector3d::UnitZ();
    Eigen::Vector3d directionOrthogonal0 = Eigen::Vector3d::UnitX();
    Eigen::Vector3d directionOrthogonal1 = Eigen::Vector3d::UnitY();
 
    Run3d(Eigen::Vector3d(6,7,3), 
          Eigen::Vector3d(6,7,4),
          Eigen::Vector3d(6,7,5), directionAligned, directionOrthogonal0, directionOrthogonal1, 2);
} 

BOOST_AUTO_TEST_CASE(TrussAlignedAngleBisector)
{
    std::cout << "*******************************************************************" << std::endl;
    std::cout << "*********** Start Example 4: Truss aligned with angle bisector ****" << std::endl;
    std::cout << "*******************************************************************" << std::endl;

    auto nodeCoord = Eigen::Vector3d(2,0,3);
    auto nodeOffset = Eigen::Vector3d(1,1,1);

    auto directionAligned = Eigen::Vector3d(1,1,1);
    auto directionOrthogonal0 = Eigen::Vector3d(-1, -1, 2);
    auto directionOrthogonal1 = Eigen::Vector3d(2, -1, -1);
 
    Run3d(nodeCoord, nodeCoord + 0.5 * nodeOffset, nodeCoord + nodeOffset, directionAligned, directionOrthogonal0, directionOrthogonal1, std::sqrt(3));
}

