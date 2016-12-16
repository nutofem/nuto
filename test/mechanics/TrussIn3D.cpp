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

#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeEnum.h"
#include "nuto/mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "nuto/mechanics/nodes/NodeEnum.h"
#include "nuto/mechanics/sections/SectionEnum.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"

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
    static const NuTo::FullVector<double, 3> mDirectionX;
    static const NuTo::FullVector<double, 3> mDirectionY;
    static const NuTo::FullVector<double, 3> mDirectionZ;
};

const NuTo::FullVector<double, 3> ParametersGeometry3D::mDirectionX = NuTo::FullVector<double, 3>::UnitX();
const NuTo::FullVector<double, 3> ParametersGeometry3D::mDirectionY = NuTo::FullVector<double, 3>::UnitY();
const NuTo::FullVector<double, 3> ParametersGeometry3D::mDirectionZ = NuTo::FullVector<double, 3>::UnitZ();

// Input:   3 node vectors (quadratic truss element),
//          1 direction vector that points in the direction of the truss (for the load)
//          2 orthogonal direction vectors (for the BC)
void Run3d(NuTo::FullVector<double, -1> rNodeCoords0, NuTo::FullVector<double, -1> rNodeCoords1, NuTo::FullVector<double, -1> rNodeCoords2, NuTo::FullVector<double, -1> rDirectionAligned, NuTo::FullVector<double, -1> rDirectionOrthogonal0, NuTo::FullVector<double, -1> rDirectionOrthogonal1)
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

    int fibreSection = myStructure.SectionCreate(NuTo::eSectionType::TRUSS);
    myStructure.SectionSetArea(fibreSection, ParametersGeometry3D::mCrossSection);

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

    myStructure.ConstraintLinearSetDisplacementNode(node0, ParametersGeometry3D::mDirectionX, 0);
    myStructure.ConstraintLinearSetDisplacementNode(node0, ParametersGeometry3D::mDirectionY, 0);
    myStructure.ConstraintLinearSetDisplacementNode(node0, ParametersGeometry3D::mDirectionZ, 0);

    myStructure.ConstraintLinearSetDisplacementNode(node1, rDirectionOrthogonal0, 0);
    myStructure.ConstraintLinearSetDisplacementNode(node2, rDirectionOrthogonal0, 0);

    myStructure.ConstraintLinearSetDisplacementNode(node1, rDirectionOrthogonal1, 0);
    myStructure.ConstraintLinearSetDisplacementNode(node2, rDirectionOrthogonal1, 0);
    //**********************************************
    //          Loads
    //**********************************************

    int load = myStructure.LoadCreateNodeForce(0, node2, rDirectionAligned, 1);

    NuTo::FullMatrix<double, 2, 2> timeDependentLoad;
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

    NuTo::FullVector<double, -1> displacements;
    boost::filesystem::remove_all(resultPath);

    myStructure.NodeGetDisplacements(node2, displacements);

    if (std::abs(displacements.norm() - std::sqrt(3)) > 1e-8)
        throw NuTo::Exception("The calculated displacements do not agree with the analytical solution!");
}

int main(int argc, char* argv[])
{
    // node coordinates
    NuTo::FullVector<double, ParametersGeometry3D::mDimension> nodeCoords0;
    NuTo::FullVector<double, ParametersGeometry3D::mDimension> nodeCoords1;
    NuTo::FullVector<double, ParametersGeometry3D::mDimension> nodeCoords2;

    // directions
    NuTo::FullVector<double, ParametersGeometry3D::mDimension> directionAligned;
    NuTo::FullVector<double, ParametersGeometry3D::mDimension> directionOrthogonal0;
    NuTo::FullVector<double, ParametersGeometry3D::mDimension> directionOrthogonal1;

    try
    {
        std::cout << "*******************************************************************" << std::endl;
        std::cout << "*********** Start Example 1: Truss aligned with x-axis ************" << std::endl;
        std::cout << "*******************************************************************" << std::endl;

        nodeCoords0[0] = 1.0;
        nodeCoords0[1] = 5.0;
        nodeCoords0[2] = 3.0;

        nodeCoords1[0] = 1.0 + 0.5 * std::sqrt(3);
        nodeCoords1[1] = 5.0;
        nodeCoords1[2] = 3.0;

        nodeCoords2[0] = 1.0 + 1.0 * std::sqrt(3);
        nodeCoords2[1] = 5.0;
        nodeCoords2[2] = 3.0;

        directionAligned = Eigen::Vector3d::UnitX();
        directionOrthogonal0 = Eigen::Vector3d::UnitY();
        directionOrthogonal1 = Eigen::Vector3d::UnitZ();

        Run3d(nodeCoords0, nodeCoords1, nodeCoords2, directionAligned, directionOrthogonal0, directionOrthogonal1);

        std::cout << "*******************************************************************" << std::endl;
        std::cout << "*********** Start Example 2: Truss aligned with y-axis ************" << std::endl;
        std::cout << "*******************************************************************" << std::endl;

        nodeCoords0[0] = 0.0;
        nodeCoords0[1] = 0.0;
        nodeCoords0[2] = 0.0;

        nodeCoords1[0] = 0.0;
        nodeCoords1[1] = 0.0 + 0.5 * std::sqrt(3);
        nodeCoords1[2] = 0.0;

        nodeCoords2[0] = 0.0;
        nodeCoords2[1] = 0.0 + 1.0 * std::sqrt(3);
        nodeCoords2[2] = 0.0;

        directionAligned = Eigen::Vector3d::UnitY();
        directionOrthogonal0 = Eigen::Vector3d::UnitX();
        directionOrthogonal1 = Eigen::Vector3d::UnitZ();

        Run3d(nodeCoords0, nodeCoords1, nodeCoords2, directionAligned, directionOrthogonal0, directionOrthogonal1);

        std::cout << "*******************************************************************" << std::endl;
        std::cout << "*********** Start Example 3: Truss aligned with z-axis ************" << std::endl;
        std::cout << "*******************************************************************" << std::endl;

        nodeCoords0[0] = 6.0;
        nodeCoords0[1] = 7.0;
        nodeCoords0[2] = 3.0;

        nodeCoords1[0] = 6.0;
        nodeCoords1[1] = 7.0;
        nodeCoords1[2] = 3.0 + 0.5 * std::sqrt(3);

        nodeCoords2[0] = 6.0;
        nodeCoords2[1] = 7.0;
        nodeCoords2[2] = 3.0 + 1.0 * std::sqrt(3);

        directionAligned = Eigen::Vector3d::UnitZ();
        directionOrthogonal0 = Eigen::Vector3d::UnitX();
        directionOrthogonal1 = Eigen::Vector3d::UnitY();

        Run3d(nodeCoords0, nodeCoords1, nodeCoords2, directionAligned, directionOrthogonal0, directionOrthogonal1);

        std::cout << "*******************************************************************" << std::endl;
        std::cout << "*********** Start Example 4: Truss aligned with angle bisector ****" << std::endl;
        std::cout << "*******************************************************************" << std::endl;

        nodeCoords0[0] = 2.0;
        nodeCoords0[1] = 0.0;
        nodeCoords0[2] = 3.0;

        nodeCoords1[0] = 2.0 + 0.5;
        nodeCoords1[1] = 0.0 + 0.5;
        nodeCoords1[2] = 3.0 + 0.5;

        nodeCoords2[0] = 2.0 + 1.0;
        nodeCoords2[1] = 0.0 + 1.0;
        nodeCoords2[2] = 3.0 + 1.0;

        directionAligned[0] = 1.0;
        directionAligned[1] = 1.0;
        directionAligned[2] = 1.0;
        directionAligned.normalize();

        directionOrthogonal0[0] = -1.0;
        directionOrthogonal0[1] = -1.0;
        directionOrthogonal0[2] = 2.0;

        directionOrthogonal1[0] = 2.0;
        directionOrthogonal1[1] = -1.0;
        directionOrthogonal1[2] = -1.0;

        Run3d(nodeCoords0, nodeCoords1, nodeCoords2, directionAligned, directionOrthogonal0, directionOrthogonal1);

    } catch (NuTo::Exception& e)
    {
        std::cout << "## Test failed ##" << std::endl;
        std::cout << e.ErrorMessage();
        return EXIT_FAILURE;
    }
    catch (...)
    {
        std::cout << "## Test failed ##" << std::endl;
        std::cout << "Something else went wrong :-(" << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "## Test successful ##" << std::endl;

}

