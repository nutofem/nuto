//============================================================================
// Name        : InterfaceElements.cpp
// Author      : Philip Huschke
// Version     : 02 Sep 2015
// Copyright   :
// Description : Test for the interface element proposed by Goodman et al.
//               Fibre pullout test:
//
//
//============================================================================
#include "BoostUnitTest.h"
#include <boost/filesystem.hpp>
#include <iostream>
#include <fstream>
#include <string>

#include "mechanics/MechanicsEnums.h"
#include "mechanics/constraints/ConstraintCompanion.h"
#include "mechanics/sections/SectionPlane.h"
#include "mechanics/sections/SectionTruss.h"
#include "mechanics/sections/SectionFibreMatrixBond.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "visualize/VisualizeEnum.h"

constexpr unsigned int dimension = 2;

class Parameters
{
public:

    static const int mDimension = dimension;

    static const bool mPerformLineSearch = true;
    static const bool mAutomaticTimeStepping = true;

    static constexpr double mMatrixYoungsModulus = 4.0e4;
    static constexpr double mMatrixPoissonsRatio = 0.2;
    static constexpr double mMatrixThickness = 0.2;

    static constexpr double mFibreYoungsModulus = 2.1e5;
    static constexpr double mFibrePoissonsRatio = 0.2;
    static constexpr double mFibreCrossSection = 0.1;
    static constexpr double mFibreCircumference = 1.1;

    static constexpr double mInterfaceNormalStiffness = 1e6;
    static constexpr double mAlpha = 1;
    static constexpr double mMaxBondStress = 4e3;
    static constexpr double mResidualBondStress = 1e3;
    static constexpr double mSlipAtMaxBondStress = 0.5;
    static constexpr double mSlipAtResidualBondStress = 5;

    static constexpr double mTimeStep = 1.0e-2;
    static constexpr double mMinTimeStep = 1.0e-5;
    static constexpr double mMaxTimeStep = 1.0e-2;
    static constexpr double mToleranceForce = 1e-6;
    static constexpr double mSimulationTime = 1.0;
    static constexpr double mLoad = 10.0;
};

//////////////////////////////////////////////////////////
//  MAIN
//////////////////////////////////////////////////////////

void Run()
{
    std::cout << "***********************************" << std::endl;
    std::cout << "**      Structure                **" << std::endl;
    std::cout << "***********************************" << std::endl;
    
    const boost::filesystem::path resultDir(boost::filesystem::initial_path().string() + "/results_interface_elements/");
    const boost::filesystem::path meshFile(boost::filesystem::initial_path().string() + "/meshes/InterfaceElements.msh");
    
    NuTo::Structure myStructure(Parameters::mDimension);
    myStructure.SetVerboseLevel(10);
    myStructure.SetShowTime(false);
    
    std::cout << "***********************************" << std::endl;
    std::cout << "**      Integration Scheme       **" << std::endl;
    std::cout << "***********************************" << std::endl;
    
    NuTo::NewmarkDirect myIntegrationScheme(&myStructure);
    myIntegrationScheme.SetTimeStep(Parameters::mTimeStep);
    myIntegrationScheme.SetMinTimeStep(Parameters::mMinTimeStep);
    myIntegrationScheme.SetMaxTimeStep(Parameters::mMaxTimeStep);
    myIntegrationScheme.SetToleranceForce(Parameters::mToleranceForce);
    myIntegrationScheme.SetAutomaticTimeStepping(Parameters::mAutomaticTimeStepping);
    myIntegrationScheme.SetPerformLineSearch(Parameters::mPerformLineSearch);
    myIntegrationScheme.SetResultDirectory(resultDir.string(), true);
    
    std::cout << "***********************************" << std::endl;
    std::cout << "**      Section                  **" << std::endl;
    std::cout << "***********************************" << std::endl;
    
    auto matrixSection = NuTo::SectionPlane::Create(Parameters::mMatrixThickness, false);
    auto fibreSection = NuTo::SectionTruss::Create(Parameters::mFibreCrossSection);
    auto fibreMatrixBond = NuTo::SectionFibreMatrixBond::Create(Parameters::mFibreCircumference);
    
    std::cout << "***********************************" << std::endl;
    std::cout << "**      Material                 **" << std::endl;
    std::cout << "***********************************" << std::endl;
    
    int matrixMaterial = myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, Parameters::mMatrixYoungsModulus);
    myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, Parameters::mMatrixPoissonsRatio);
    
    int fibreMaterial = myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    myStructure.ConstitutiveLawSetParameterDouble(fibreMaterial, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, Parameters::mFibreYoungsModulus);
    myStructure.ConstitutiveLawSetParameterDouble(fibreMaterial, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, Parameters::mFibrePoissonsRatio);
    
    int interfaceMaterial = myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::FIBRE_MATRIX_BOND_STRESS_SLIP);
    myStructure.ConstitutiveLawSetParameterDouble(interfaceMaterial, NuTo::Constitutive::eConstitutiveParameter::NORMAL_STIFFNESS, Parameters::mInterfaceNormalStiffness);
    myStructure.ConstitutiveLawSetParameterDouble(interfaceMaterial, NuTo::Constitutive::eConstitutiveParameter::ALPHA, Parameters::mAlpha);
    myStructure.ConstitutiveLawSetParameterDouble(interfaceMaterial, NuTo::Constitutive::eConstitutiveParameter::MAX_BOND_STRESS, Parameters::mMaxBondStress);
    myStructure.ConstitutiveLawSetParameterDouble(interfaceMaterial, NuTo::Constitutive::eConstitutiveParameter::RESIDUAL_BOND_STRESS, Parameters::mResidualBondStress);
    myStructure.ConstitutiveLawSetParameterDouble(interfaceMaterial, NuTo::Constitutive::eConstitutiveParameter::SLIP_AT_MAX_BOND_STRESS, Parameters::mSlipAtMaxBondStress);
    myStructure.ConstitutiveLawSetParameterDouble(interfaceMaterial, NuTo::Constitutive::eConstitutiveParameter::SLIP_AT_RESIDUAL_BOND_STRESS, Parameters::mSlipAtResidualBondStress);
    
    std::cout << "***********************************" << std::endl;
    std::cout << "**      Interpolation Type       **" << std::endl;
    std::cout << "***********************************" << std::endl;
    
    int matrixInterpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TRIANGLE2D);
    myStructure.InterpolationTypeAdd(matrixInterpolationType, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    myStructure.InterpolationTypeAdd(matrixInterpolationType, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    
    int fibreInterpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TRUSSXD);
    myStructure.InterpolationTypeAdd(fibreInterpolationType, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    myStructure.InterpolationTypeAdd(fibreInterpolationType, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    
    int interfaceInterpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::INTERFACE);
    myStructure.InterpolationTypeAdd(interfaceInterpolationType, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    myStructure.InterpolationTypeAdd(interfaceInterpolationType, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    
    std::cout << "***********************************" << std::endl;
    std::cout << "**      Import Mesh File         **" << std::endl;
    std::cout << "***********************************" << std::endl;
    
    auto createdGroupIdMatrix = myStructure.ImportFromGmsh(meshFile.string());
    
    int groupIdFibre = createdGroupIdMatrix[0].first;
    int groupIdMatrix = createdGroupIdMatrix[1].first;
    
    std::cout << "***********************************" << std::endl;
    std::cout << "**      Matrix                   **" << std::endl;
    std::cout << "***********************************" << std::endl;
    
    myStructure.ElementGroupSetInterpolationType(groupIdMatrix, matrixInterpolationType);
    myStructure.ElementGroupSetSection(groupIdMatrix, matrixSection);
    myStructure.ElementGroupSetConstitutiveLaw(groupIdMatrix, matrixMaterial);
    
    std::cout << "***********************************" << std::endl;
    std::cout << "**      Fibre                    **" << std::endl;
    std::cout << "***********************************" << std::endl;
    
    myStructure.ElementGroupSetInterpolationType(groupIdFibre, fibreInterpolationType);
    myStructure.ElementGroupSetSection(groupIdFibre, fibreSection);
    myStructure.ElementGroupSetConstitutiveLaw(groupIdFibre, fibreMaterial);
    
    myStructure.ElementTotalConvertToInterpolationType(1e-6, 10);
    
    std::cout << "***********************************" << std::endl;
    std::cout << "**      Interface                **" << std::endl;
    std::cout << "***********************************" << std::endl;
    
    auto pairGroupFiberGroupBond = myStructure.InterfaceElementsCreate(groupIdFibre, interfaceInterpolationType, fibreInterpolationType);
    
    int groupEleFiber   = pairGroupFiberGroupBond.first;
    int groupEleBond    = pairGroupFiberGroupBond.second;
    
    myStructure.ElementGroupSetConstitutiveLaw(groupEleFiber, fibreMaterial);
    myStructure.ElementGroupSetSection(groupEleFiber, fibreSection);
    
    myStructure.ElementGroupSetConstitutiveLaw(groupEleBond, interfaceMaterial);
    myStructure.ElementGroupSetSection(groupEleBond, fibreMatrixBond);
    
    std::cout << "***********************************" << std::endl;
    std::cout << "**      Boundary Conditions      **" << std::endl;
    std::cout << "***********************************" << std::endl;
    
    
    auto& groupNodeBCLeft = myStructure.GroupGetNodesAtCoordinate(NuTo::eDirection::X, 0);
    myStructure.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS,
            NuTo::Constraint::Component(groupNodeBCLeft, {NuTo::eDirection::X}));
    
    auto& nodeLeft = myStructure.NodeGetAtCoordinate(Eigen::Vector2d(0, 5));
    myStructure.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS,
            NuTo::Constraint::Component(nodeLeft, {NuTo::eDirection::Y}));
    
    auto& nodeFibre = *myStructure.NodeGetNodePtr(205); // 205 is the node id of the fibre
    auto loadFunction = [=](double time)
    {
        if (time < 0.5 * Parameters::mSimulationTime)
            return 2 * time * Parameters::mLoad / Parameters::mSimulationTime;
        return Parameters::mLoad * (2. - 2 * time / Parameters::mSimulationTime);
    };
    
    myStructure.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS,
            NuTo::Constraint::Component(nodeFibre, {NuTo::eDirection::Y}));
    
    myStructure.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS,
            NuTo::Constraint::Component(nodeFibre, {NuTo::eDirection::X}, loadFunction));
    
    std::cout << "***********************************" << std::endl;
    std::cout << "**      Visualization            **" << std::endl;
    std::cout << "***********************************" << std::endl;
    
    int visualizationGroup = myStructure.GroupCreate(NuTo::eGroupId::Elements);
    myStructure.GroupAddElementsTotal(visualizationGroup);
    
    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::DISPLACEMENTS);
    
    
    std::cout << "***********************************" << std::endl;
    std::cout << "**      Solver                   **" << std::endl;
    std::cout << "***********************************" << std::endl;
    
    myStructure.NodeBuildGlobalDofs();
    myStructure.CalculateMaximumIndependentSets();
    
    myIntegrationScheme.Solve(Parameters::mSimulationTime);
    
    std::cout << "Results written to " + resultDir.string() << std::endl;
}

BOOST_AUTO_TEST_CASE(InterfaceElements)
{
    BOOST_CHECK_NO_THROW(Run());
}
