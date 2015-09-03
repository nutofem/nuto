//============================================================================
// Name        : InterfaceElements.cpp
// Author      : Philip Huschke
// Version     : 02 Sep 2015
// Copyright   :
// Description : Test for the interface element proposed by Goodman et al.
//               Fibre pullout test:
//
//                    x-----x-----x
//                    |     |     |
//                    |  3  |  4  |           matrix element
//                    |     |     |
//                    x-----x-----x
//                    |     |     |
//                    |  7  |  8  |           interface element
//                    |     |     |
//                    x==9==x==10=x  -> F     fibre element
//                    |     |     |
//                    |  5  |  6  |           interface element
//                    |     |     |
//                    x-----x-----x
//                    |     |     |
//                    |  1  |  2  |           matrix element
//                    |     |     |
//                    x-----x-----x
//
//============================================================================



#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"

#include <boost-1_55/boost/filesystem.hpp>
#include <boost-1_55/boost/lexical_cast.hpp>
#include <boost-1_55/boost/progress.hpp>

#include <iostream>
#include <fstream>
#include <string>

const unsigned int dimension = 2;

class Parameters
{
public:

    static constexpr int mDimension = dimension;

    static constexpr double mMatrixYoungsModulus = 2.0e3;
    static constexpr double mMatrixPoissonsRatio = 0.3;
    static constexpr double mMatrixThickness = 1.0;

    static constexpr double mInterfaceYoungsModulus = 1.0e3;

    static constexpr double mInterfaceNormalStiffness = 1e10;
    static constexpr double mInterfaceTangentialStiffness = 1e4;

    static constexpr double mFibreYoungsModulus = 2.0e4;
    static constexpr double mFibreCrossSection = 1.0;

    static constexpr double mTimeStep = 1.e0;
    static constexpr double mToleranceForce = 1e-6;
    static constexpr double mSimulationTime = 1.0;
    static constexpr double mLoad = 1.0;
    static constexpr bool mAutomaticTimeStepping = false;

    static const NuTo::FullVector<double, dimension> mDirectionX;
    static const NuTo::FullVector<double, dimension> mDirectionY;
};

const NuTo::FullVector<double, dimension> Parameters::mDirectionX = NuTo::FullVector<double, 2>::UnitX();
const NuTo::FullVector<double, dimension> Parameters::mDirectionY = NuTo::FullVector<double, 2>::UnitY();

int main()
{

    std::cout << "**********************************************" << std::endl;
    std::cout << "                  Structure                   " << std::endl;
    std::cout << "**********************************************" << std::endl;

    NuTo::Structure myStructure(Parameters::mDimension);
    myStructure.SetVerboseLevel(10);
    myStructure.SetShowTime(false);

    std::cout << "**********************************************" << std::endl;
    std::cout << "                  Integration Scheme          " << std::endl;
    std::cout << "**********************************************" << std::endl;

    boost::filesystem::path resultPath(boost::filesystem::initial_path().string() + "/result_test_interface_elements/");

    NuTo::NewmarkDirect myIntegrationScheme(&myStructure);
    myIntegrationScheme.SetTimeStep(Parameters::mTimeStep);
    myIntegrationScheme.SetToleranceForce(Parameters::mToleranceForce);
    myIntegrationScheme.SetAutomaticTimeStepping(Parameters::mAutomaticTimeStepping);
    myIntegrationScheme.SetResultDirectory(resultPath.string(), true);

    std::cout << "**********************************************" << std::endl;
    std::cout << "                  Section                     " << std::endl;
    std::cout << "**********************************************" << std::endl;

    int interfaceSection = myStructure.SectionCreate(NuTo::Section::PLANE_STRAIN);

    int matrixSection = myStructure.SectionCreate(NuTo::Section::PLANE_STRAIN);
    myStructure.SectionSetThickness(matrixSection, Parameters::mMatrixThickness);

    int fibreSection = myStructure.SectionCreate(NuTo::Section::TRUSS);
    myStructure.SectionSetArea(fibreSection, Parameters::mFibreCrossSection);

    std::cout << "**********************************************" << std::endl;
    std::cout << "                  Material                    " << std::endl;
    std::cout << "**********************************************" << std::endl;

    int fibreMaterial = myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    myStructure.ConstitutiveLawSetParameterDouble(fibreMaterial, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, Parameters::mFibreYoungsModulus);

    int matrixMaterial = myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, Parameters::mMatrixYoungsModulus);
    myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, Parameters::mMatrixPoissonsRatio);

    int interfaceMaterial = myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::INTERFACE_GOODMAN);
    myStructure.ConstitutiveLawSetParameterDouble(interfaceMaterial, NuTo::Constitutive::eConstitutiveParameter::NORMAL_STIFFNESS, Parameters::mInterfaceNormalStiffness);
    myStructure.ConstitutiveLawSetParameterDouble(interfaceMaterial, NuTo::Constitutive::eConstitutiveParameter::TANGENTIAL_STIFFNESS, Parameters::mInterfaceTangentialStiffness);

    std::cout << "**********************************************" << std::endl;
    std::cout << "                  Interpolation Type          " << std::endl;
    std::cout << "**********************************************" << std::endl;

    int fibreInterpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TRUSSXD);
    myStructure.InterpolationTypeAdd(fibreInterpolationType, NuTo::Node::COORDINATES, NuTo::Interpolation::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(fibreInterpolationType, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::EQUIDISTANT1);

    int matrixInterpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::QUAD2D);
    myStructure.InterpolationTypeAdd(matrixInterpolationType, NuTo::Node::COORDINATES, NuTo::Interpolation::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(matrixInterpolationType, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::EQUIDISTANT1);

    int interfaceInterpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::INTERFACE);
    myStructure.InterpolationTypeAdd(interfaceInterpolationType, NuTo::Node::COORDINATES, NuTo::Interpolation::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(interfaceInterpolationType, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::EQUIDISTANT1);

    std::cout << "**********************************************" << std::endl;
    std::cout << "                  Nodes                       " << std::endl;
    std::cout << "**********************************************" << std::endl;

    NuTo::FullVector<double, Eigen::Dynamic> nodeCoords(Parameters::mDimension);
    std::set<NuTo::Node::eAttributes> dofs;
    dofs.insert(NuTo::Node::COORDINATES);
    dofs.insert(NuTo::Node::DISPLACEMENTS);

    int node01 = myStructure.NodeCreate(Eigen::Vector2d(0.0, 0.0), dofs);
    int node02 = myStructure.NodeCreate(Eigen::Vector2d(1.0, 0.0), dofs);
    int node03 = myStructure.NodeCreate(Eigen::Vector2d(2.0, 0.0), dofs);

    int node04 = myStructure.NodeCreate(Eigen::Vector2d(0.0, 1.0), dofs);
    int node05 = myStructure.NodeCreate(Eigen::Vector2d(1.0, 1.0), dofs);
    int node06 = myStructure.NodeCreate(Eigen::Vector2d(2.0, 1.0), dofs);

    int node07 = myStructure.NodeCreate(Eigen::Vector2d(0.0, 1.0), dofs);
    int node08 = myStructure.NodeCreate(Eigen::Vector2d(1.0, 1.0), dofs);
    int node09 = myStructure.NodeCreate(Eigen::Vector2d(2.0, 1.0), dofs);

    int node10 = myStructure.NodeCreate(Eigen::Vector2d(0.0, 1.0), dofs);
    int node11 = myStructure.NodeCreate(Eigen::Vector2d(1.0, 1.0), dofs);
    int node12 = myStructure.NodeCreate(Eigen::Vector2d(2.0, 1.0), dofs);

    int node13 = myStructure.NodeCreate(Eigen::Vector2d(0.0, 2.0), dofs);
    int node14 = myStructure.NodeCreate(Eigen::Vector2d(1.0, 2.0), dofs);
    int node15 = myStructure.NodeCreate(Eigen::Vector2d(2.0, 2.0), dofs);

    std::cout << "**********************************************" << std::endl;
    std::cout << "                  Matrix Geometry             " << std::endl;
    std::cout << "**********************************************" << std::endl;

    NuTo::FullVector<int, Eigen::Dynamic> nodeIndicesMatrix(4);

    nodeIndicesMatrix[0] = node01;
    nodeIndicesMatrix[1] = node02;
    nodeIndicesMatrix[2] = node05;
    nodeIndicesMatrix[3] = node04;
    int elementMatirx01 = myStructure.ElementCreate(matrixInterpolationType, nodeIndicesMatrix, NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::NOIPDATA);

    nodeIndicesMatrix[0] = node02;
    nodeIndicesMatrix[1] = node03;
    nodeIndicesMatrix[2] = node06;
    nodeIndicesMatrix[3] = node05;
    int elementMatirx02 = myStructure.ElementCreate(matrixInterpolationType, nodeIndicesMatrix, NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::NOIPDATA);

    nodeIndicesMatrix[0] = node10;
    nodeIndicesMatrix[1] = node11;
    nodeIndicesMatrix[2] = node14;
    nodeIndicesMatrix[3] = node13;
    int elementMatirx03 = myStructure.ElementCreate(matrixInterpolationType, nodeIndicesMatrix, NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::NOIPDATA);

    nodeIndicesMatrix[0] = node11;
    nodeIndicesMatrix[1] = node12;
    nodeIndicesMatrix[2] = node15;
    nodeIndicesMatrix[3] = node14;
    int elementMatirx04 = myStructure.ElementCreate(matrixInterpolationType, nodeIndicesMatrix, NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::NOIPDATA);

    int groupElementsMatrix = myStructure.GroupCreate(NuTo::Groups::eGroupId::Elements);
    myStructure.GroupAddElement(groupElementsMatrix, elementMatirx01);
    myStructure.GroupAddElement(groupElementsMatrix, elementMatirx02);
    myStructure.GroupAddElement(groupElementsMatrix, elementMatirx03);
    myStructure.GroupAddElement(groupElementsMatrix, elementMatirx04);

    myStructure.ElementGroupSetSection(groupElementsMatrix, matrixMaterial);
    myStructure.ElementGroupSetConstitutiveLaw(groupElementsMatrix, matrixSection);

    std::cout << "**********************************************" << std::endl;
    std::cout << "                  Interface Geometry          " << std::endl;
    std::cout << "**********************************************" << std::endl;

    NuTo::FullVector<int, Eigen::Dynamic> nodeIndicesInterface(4);

    nodeIndicesInterface[0] = node04;
    nodeIndicesInterface[1] = node05;
    nodeIndicesInterface[2] = node08;
    nodeIndicesInterface[3] = node07;
    int elementInterface05 = myStructure.ElementCreate(interfaceInterpolationType, nodeIndicesInterface, NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::NOIPDATA);

    nodeIndicesInterface[0] = node05;
    nodeIndicesInterface[1] = node06;
    nodeIndicesInterface[2] = node09;
    nodeIndicesInterface[3] = node08;
    int elementInterface06 = myStructure.ElementCreate(interfaceInterpolationType, nodeIndicesInterface, NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::NOIPDATA);

    nodeIndicesInterface[0] = node07;
    nodeIndicesInterface[1] = node08;
    nodeIndicesInterface[2] = node11;
    nodeIndicesInterface[3] = node10;
    int elementInterface07 = myStructure.ElementCreate(interfaceInterpolationType, nodeIndicesInterface, NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::NOIPDATA);

    nodeIndicesInterface[0] = node08;
    nodeIndicesInterface[1] = node09;
    nodeIndicesInterface[2] = node12;
    nodeIndicesInterface[3] = node11;
    int elementInterface08 = myStructure.ElementCreate(interfaceInterpolationType, nodeIndicesInterface, NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::NOIPDATA);

    int groupElementsInterface = myStructure.GroupCreate(NuTo::Groups::eGroupId::Elements);
    myStructure.GroupAddElement(groupElementsInterface, elementInterface05);
    myStructure.GroupAddElement(groupElementsInterface, elementInterface06);
    myStructure.GroupAddElement(groupElementsInterface, elementInterface07);
    myStructure.GroupAddElement(groupElementsInterface, elementInterface08);

    myStructure.ElementGroupSetSection(groupElementsInterface, interfaceSection);
    myStructure.ElementGroupSetConstitutiveLaw(groupElementsInterface, interfaceMaterial);

    myStructure.NodeBuildGlobalDofs();

    std::cout << "**********************************************" << std::endl;
    std::cout << "                  Fibre Geometry              " << std::endl;
    std::cout << "**********************************************" << std::endl;

    NuTo::FullVector<int, Eigen::Dynamic> nodeIndicesFibre(2);

    nodeIndicesFibre[0] = node07;
    nodeIndicesFibre[1] = node08;
    int elementFibre09 = myStructure.ElementCreate(fibreInterpolationType, nodeIndicesFibre, NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::NOIPDATA);

    nodeIndicesFibre[0] = node08;
    nodeIndicesFibre[1] = node09;
    int elementFibre10 = myStructure.ElementCreate(fibreInterpolationType, nodeIndicesFibre, NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::NOIPDATA);

    int groupElementsFibre = myStructure.GroupCreate(NuTo::Groups::eGroupId::Elements);
    myStructure.GroupAddElement(groupElementsFibre, elementFibre09);
    myStructure.GroupAddElement(groupElementsFibre, elementFibre10);

    myStructure.ElementGroupSetSection(groupElementsFibre, fibreSection);
    myStructure.ElementGroupSetConstitutiveLaw(groupElementsFibre, fibreMaterial);

    myStructure.NodeBuildGlobalDofs();

    std::cout << "**********************************************" << std::endl;
    std::cout << "                  Boundary Conditions         " << std::endl;
    std::cout << "**********************************************" << std::endl;

    int groupNodeBCBottom = myStructure.GroupCreate(NuTo::Groups::eGroupId::Nodes);
    myStructure.GroupAddNode(groupNodeBCBottom, node01);
    myStructure.GroupAddNode(groupNodeBCBottom, node02);
    myStructure.GroupAddNode(groupNodeBCBottom, node03);
    myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodeBCBottom, Parameters::mDirectionY, 0);

    int groupNodeBCLeft = myStructure.GroupCreate(NuTo::Groups::eGroupId::Nodes);
    myStructure.GroupAddNode(groupNodeBCLeft, node01);
    myStructure.GroupAddNode(groupNodeBCLeft, node04);
    myStructure.GroupAddNode(groupNodeBCLeft, node10);
    myStructure.GroupAddNode(groupNodeBCLeft, node13);
    myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodeBCLeft, Parameters::mDirectionX, 0);

    int groupNodeBCRight = myStructure.GroupCreate(NuTo::Groups::eGroupId::Nodes);
    myStructure.GroupAddNode(groupNodeBCRight, node03);
    myStructure.GroupAddNode(groupNodeBCRight, node06);
    myStructure.GroupAddNode(groupNodeBCRight, node12);
    myStructure.GroupAddNode(groupNodeBCRight, node15);
    myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodeBCRight, Parameters::mDirectionX, 0);

    myStructure.NodeBuildGlobalDofs();

    std::cout << "**********************************************" << std::endl;
    std::cout << "                  Loads                       " << std::endl;
    std::cout << "**********************************************" << std::endl;

    int loadCase = 0;
    myStructure.LoadCreateNodeForce(loadCase, node09, Parameters::mDirectionX, 12e0);

    myStructure.LoadCreateNodeForce(loadCase, node13, Parameters::mDirectionY, -5e0);
    myStructure.LoadCreateNodeForce(loadCase, node14, Parameters::mDirectionY, -10e0);
    myStructure.LoadCreateNodeForce(loadCase, node15, Parameters::mDirectionY, -5e0);

    std::cout << "**********************************************" << std::endl;
    std::cout << "                  Visualisation               " << std::endl;
    std::cout << "**********************************************" << std::endl;

    myStructure.AddVisualizationComponentDisplacements();
    myStructure.AddVisualizationComponentSection();
    myStructure.AddVisualizationComponentConstitutive();
    myStructure.AddVisualizationComponentElement();

    std::cout << "**********************************************" << std::endl;
    std::cout << "                  Solver                      " << std::endl;
    std::cout << "**********************************************" << std::endl;

    myStructure.NodeBuildGlobalDofs();
    myStructure.CalculateMaximumIndependentSets();

    NuTo::FullMatrix<double, 2, 2> timeDependentLoad;
    timeDependentLoad(0, 0) = 0;
    timeDependentLoad(1, 0) = Parameters::mSimulationTime;
    timeDependentLoad(0, 1) = 0;
    timeDependentLoad(1, 1) = Parameters::mLoad;

    myIntegrationScheme.SetTimeDependentLoadCase(loadCase, timeDependentLoad);
    myIntegrationScheme.Solve(Parameters::mSimulationTime);

    std::cout << "**********************************************" << std::endl;
    std::cout << "                  End                         " << std::endl;
    std::cout << "**********************************************" << std::endl;
}

