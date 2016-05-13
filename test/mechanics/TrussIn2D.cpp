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

#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"

#include <boost/filesystem.hpp>
#include <iostream>
#include <fstream>
#include <string>

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
    static constexpr int mDimension = 2;
    static constexpr double mCrossSection = 0.1;
    static const NuTo::FullVector<double, 2> mDirectionX;
    static const NuTo::FullVector<double, 2> mDirectionY;
};

const NuTo::FullVector<double, 2> ParametersGeometry2D::mDirectionX = NuTo::FullVector<double, 2>::UnitX();
const NuTo::FullVector<double, 2> ParametersGeometry2D::mDirectionY = NuTo::FullVector<double, 2>::UnitY();

void Run2d(NuTo::FullVector<double, -1> rNodeCoords0, NuTo::FullVector<double, -1> rNodeCoords1, NuTo::FullVector<double, -1> rNodeCoords2, NuTo::FullVector<double, -1> rDirectionAligned, NuTo::FullVector<double, -1> rDirectionOrthogonal)
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
    myIntegrationScheme.SetResultDirectory(resultPath.string(), false);

    //**********************************************
    //          Section
    //**********************************************

    int fibreSection = myStructure.SectionCreate(NuTo::Section::TRUSS);
    myStructure.SectionSetArea(fibreSection, ParametersGeometry2D::mCrossSection);

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
    myStructure.InterpolationTypeAdd(fibreInterpolationType, NuTo::Node::COORDINATES, NuTo::Interpolation::EQUIDISTANT2);
    myStructure.InterpolationTypeAdd(fibreInterpolationType, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::EQUIDISTANT2);
    myStructure.InterpolationTypeSetIntegrationType(fibreInterpolationType, NuTo::IntegrationType::IntegrationType1D2NGauss2Ip, NuTo::IpData::STATICDATA);

    //**********************************************
    //          Geometry
    //**********************************************

    // Nodes
    int node0 = myStructure.NodeCreate(rNodeCoords0);
    int node1 = myStructure.NodeCreate(rNodeCoords1);
    int node2 = myStructure.NodeCreate(rNodeCoords2);

    NuTo::FullVector<int, -1> nodeIndices(3);
    nodeIndices[0] = node0;
    nodeIndices[1] = node1;
    nodeIndices[2] = node2;
    myStructure.ElementCreate(fibreInterpolationType, nodeIndices, NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::NOIPDATA);

    myStructure.ElementTotalConvertToInterpolationType(1e-6, 10);
    myStructure.ElementTotalSetSection(fibreSection);
    myStructure.ElementTotalSetConstitutiveLaw(fibreMaterial);

    //**********************************************
    //          Boundary Conditions
    //**********************************************
    myStructure.ConstraintLinearSetDisplacementNode(node0, ParametersGeometry2D::mDirectionX, 0);
    myStructure.ConstraintLinearSetDisplacementNode(node0, ParametersGeometry2D::mDirectionY, 0);

    myStructure.ConstraintLinearSetDisplacementNode(node1, rDirectionOrthogonal, 0);
    myStructure.ConstraintLinearSetDisplacementNode(node2, rDirectionOrthogonal, 0);

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

    if (std::abs(displacements.norm() - std::sqrt(2)) > 1e-8)
        throw NuTo::Exception("The calculated displacements do not agree with the analytical solution!");
}

int main(int argc, char* argv[])
{

    // node coordinates
    NuTo::FullVector<double, -1> nodeCoords0(ParametersGeometry2D::mDimension);
    NuTo::FullVector<double, -1> nodeCoords1(ParametersGeometry2D::mDimension);
    NuTo::FullVector<double, -1> nodeCoords2(ParametersGeometry2D::mDimension);

    // directions
    NuTo::FullVector<double, -1> directionAligned(ParametersGeometry2D::mDimension);
    NuTo::FullVector<double, -1> directionOrthogonal(ParametersGeometry2D::mDimension);

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

    } catch (NuTo::Exception& e)
    {
        std::cout << "## Test failed ##" << std::endl;
        std::cout << e.ErrorMessage();
        return EXIT_FAILURE;
    }

    std::cout << "## Test successful ##" << std::endl;

}

