//============================================================================
// Name        : ConstraintsNodeToElement.cpp
// Author      : Philip Huschke
// Version     : 11 Jan 2016
// Copyright   :
// Description : Test for the constraints between nodes inside an element
//
//
//============================================================================
#include "BoostUnitTest.h"

#include "math/MathException.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/constitutive/damageLaws/DamageLawExponential.h"
#include "mechanics/elements/ElementBase.h"
#include "mechanics/integrationtypes/IntegrationTypeEnum.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "mechanics/groups/GroupEnum.h"
#include "mechanics/groups/Group.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/sections/SectionPlane.h"
#include "mechanics/sections/SectionTruss.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "visualize/VisualizeEnum.h"
#include "mechanics/constraints/Constraints.h"
#include "mechanics/constraints/ConstraintCompanion.h"

#include <boost/filesystem.hpp>
#include <iostream>
#include <fstream>
#include <string>



class Parameters
{
public:

    static constexpr bool mPerformLineSearch = true;
    static constexpr bool mAutomaticTimeStepping = true;

    static constexpr double mMatrixYoungsModulus = 4.0e4;   // concrete
    static constexpr double mMatrixPoissonsRatio = 0.2;
    static constexpr double mMatrixThickness = 2;
    static constexpr double mMatrixNonlocalRadius = 2;
    static constexpr double mMatrixTensileStrength = 3;
    static constexpr double mMatrixCompressiveStrength = 30;
    static constexpr double mMatrixFractureEnergy = 1;

    static constexpr double mFibreYoungsModulus = 2.1e5;
    static constexpr double mFibrePoissonsRatio = 0.2;
    static constexpr double mFibreCrossSection = 0.1;

    static constexpr double mTimeStep = 1.0e-1;
    static constexpr double mMinTimeStep = 1.0e-5;
    static constexpr double mMaxTimeStep = 1.0e-0;
    static constexpr double mToleranceForce = 1e-6;
    static constexpr double mSimulationTime = 1.0;
    static constexpr double mLoad = 0.01;

};


//////////////////////////////////////////////////////////
//  RUN 2D
//////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(run2d)
{
        constexpr int dimension = 2;

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Structure                **" << std::endl;
        std::cout << "***********************************" << std::endl;

        const boost::filesystem::path resultDir             (boost::filesystem::initial_path().string() + "/results_ConstraintsNodeToElement2d/");
        const boost::filesystem::path meshFilePathMatrix    (boost::filesystem::initial_path().string() + "/ConstraintsNodeToElementMatrixMesh2d.msh");
        const boost::filesystem::path meshFilePathFiber     (boost::filesystem::initial_path().string() + "/ConstraintsNodeToElementFiberMesh2d.msh");

        NuTo::Structure s(dimension);
        s.SetVerboseLevel(10);
        s.SetShowTime(false);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Integration Scheme       **" << std::endl;
        std::cout << "***********************************" << std::endl;

        NuTo::NewmarkDirect myIntegrationScheme(&s);
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

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Material                 **" << std::endl;
        std::cout << "***********************************" << std::endl;

        int matrixMaterial = s.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
        s.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, Parameters::mMatrixYoungsModulus);
        s.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, Parameters::mMatrixPoissonsRatio);

//        int matrixMaterial = s.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::GRADIENT_DAMAGE_ENGINEERING_STRESS);
//        s.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, Parameters::mMatrixYoungsModulus);
//        s.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, Parameters::mMatrixPoissonsRatio);
//        s.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::TENSILE_STRENGTH, Parameters::mMatrixTensileStrength);
//        s.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::COMPRESSIVE_STRENGTH, Parameters::mMatrixCompressiveStrength);
//        s.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS, Parameters::mMatrixNonlocalRadius);
//        s.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::FRACTURE_ENERGY, Parameters::mMatrixFractureEnergy);

        int fibreMaterial = s.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
        s.ConstitutiveLawSetParameterDouble(fibreMaterial, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, Parameters::mFibreYoungsModulus);
        s.ConstitutiveLawSetParameterDouble(fibreMaterial, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, Parameters::mFibrePoissonsRatio);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Interpolation Type       **" << std::endl;
        std::cout << "***********************************" << std::endl;

        int matrixInterpolationType = s.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TRIANGLE2D);
        s.InterpolationTypeAdd(matrixInterpolationType, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
        s.InterpolationTypeAdd(matrixInterpolationType, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
//        s.InterpolationTypeAdd(matrixInterpolationType, NuTo::Node::eDof::NONLOCALEQSTRAIN, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);


        int fibreInterpolationType = s.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TRUSSXD);
        s.InterpolationTypeAdd(fibreInterpolationType, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
        s.InterpolationTypeAdd(fibreInterpolationType, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Import Matrix Mesh       **" << std::endl;
        std::cout << "***********************************" << std::endl;

        auto createdGroupIdMatrix = s.ImportFromGmsh(meshFilePathMatrix.string());
        int groupIdMatrix = createdGroupIdMatrix[0].first;


        s.ElementGroupSetInterpolationType(groupIdMatrix, matrixInterpolationType);
        s.InterpolationTypeSetIntegrationType(matrixInterpolationType, NuTo::eIntegrationType::IntegrationType2D3NGauss3Ip);
        s.ElementTotalConvertToInterpolationType();
        s.ElementGroupSetSection(groupIdMatrix, matrixSection);
        s.ElementGroupSetConstitutiveLaw(groupIdMatrix, matrixMaterial);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Import Fiber Mesh        **" << std::endl;
        std::cout << "***********************************" << std::endl;

        auto createdGroupIdFiber = s.ImportFromGmsh(meshFilePathFiber.string());
        int groupIdFiber = createdGroupIdFiber[0].first;

        s.ElementGroupSetInterpolationType(groupIdFiber, fibreInterpolationType);
        s.InterpolationTypeSetIntegrationType(fibreInterpolationType, NuTo::eIntegrationType::IntegrationType1D2NGauss3Ip);
        s.ElementConvertToInterpolationType(groupIdFiber);
        s.ElementGroupSetSection(groupIdFiber, fibreSection);
        s.ElementGroupSetConstitutiveLaw(groupIdFiber, fibreMaterial);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Boundary Conditions      **" << std::endl;
        std::cout << "***********************************" << std::endl;

        int groupNodeBCLeft = s.GroupCreate(NuTo::eGroupId::Nodes);
        s.GroupAddNodeCoordinateRange(groupNodeBCLeft, 0, 0.0 - 1e-6, 0.0 + 1e-6);
        auto nodeGroup = s.GroupGetGroupPtr(groupNodeBCLeft)->AsGroupNode();
        s.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Constraint::Component(*nodeGroup, {NuTo::eDirection::X}));

        auto nodeLeft = s.NodeGetNodePtr(s.NodeGetIdAtCoordinate(Eigen::Vector2d(0,0), 1e-6));
        s.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Constraint::Component(*nodeLeft, {NuTo::eDirection::Y}));

        auto nodeRight = s.NodeGetNodePtr(s.NodeGetIdAtCoordinate(Eigen::Vector2d(10,0), 1e-6));
        s.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Constraint::Component(*nodeRight, {NuTo::eDirection::Y}));


        std::cout << "***********************************" << std::endl;
        std::cout << "**      Constraints              **" << std::endl;
        std::cout << "***********************************" << std::endl;


        int groupNodesFiber = s.GroupCreate(NuTo::eGroupId::Nodes);
        s.GroupAddNodesFromElements(groupNodesFiber, groupIdFiber);


        for (int nodeId : s.GroupGetMemberIds(groupNodesFiber))
            s.ConstraintLinearEquationNodeToElementCreate(nodeId, groupIdMatrix, NuTo::Node::eDof::DISPLACEMENTS);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Loads                    **" << std::endl;
        std::cout << "***********************************" << std::endl;

        int groupNodeLoadRight = s.GroupCreate(NuTo::eGroupId::Nodes);
        s.GroupAddNodeCoordinateRange(groupNodeLoadRight, 0, 10.0 - 1e-6, 10.0 + 1e-6);
        auto nodeGroupBC = s.GroupGetGroupPtr(groupNodeLoadRight)->AsGroupNode();

        s.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Constraint::Component(*nodeGroupBC, {NuTo::eDirection::X}, NuTo::Constraint::RhsRamp(Parameters::mSimulationTime, Parameters::mLoad)));

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Visualization            **" << std::endl;
        std::cout << "***********************************" << std::endl;

        s.AddVisualizationComponent(groupIdMatrix, NuTo::eVisualizeWhat::DISPLACEMENTS);
        s.AddVisualizationComponent(groupIdFiber, NuTo::eVisualizeWhat::DISPLACEMENTS);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Solver                   **" << std::endl;
        std::cout << "***********************************" << std::endl;

        s.NodeBuildGlobalDofs();
        s.CalculateMaximumIndependentSets();


        myIntegrationScheme.Solve(Parameters::mSimulationTime);



        std::cout << "***********************************" << std::endl;
        std::cout << "**      Postprocessing           **" << std::endl;
        std::cout << "***********************************" << std::endl;



        // Calculate the displacements at one of the constrained points in the matrix and the fiber and check if they match
        auto elementPtrMatrix = s.ElementGetElementPtr(0);

        // natural coordinates of (2.5, 0.25) in tetrahedron 0

        auto dispInMatrix = elementPtrMatrix->InterpolateDofGlobal(Eigen::Vector2d(0.5, 0.25), NuTo::Node::eDof::DISPLACEMENTS);
        auto coordsInMatrix = elementPtrMatrix->InterpolateDofGlobal(Eigen::Vector2d(0.5, 0.25), NuTo::Node::eDof::COORDINATES);
        std::cout << "dispInMatrix \n" << dispInMatrix << std::endl;
        std::cout << "coordsInMatrix \n" << coordsInMatrix << std::endl;

        auto elementPtrFiber = s.ElementGetElementPtr(4);

        // natural coordinates of (2.5, 0.25) in fiber 4
        Eigen::VectorXd fibreNaturalCoordiante = Eigen::VectorXd::Constant(1, -1);
        auto dispInFiber = elementPtrFiber->InterpolateDofGlobal(fibreNaturalCoordiante, NuTo::Node::eDof::DISPLACEMENTS);
        auto coordsInFiber = elementPtrFiber->InterpolateDofGlobal(fibreNaturalCoordiante, NuTo::Node::eDof::COORDINATES);
        std::cout << "dispInFiber \n" << dispInFiber << std::endl;
        std::cout << "coordsInFiber \n" << coordsInFiber << std::endl;

        BoostUnitTest::CheckVector(dispInFiber, dispInMatrix, 2, 1.e-6);
        BoostUnitTest::CheckVector(coordsInFiber, coordsInMatrix, 2, 1.e-6);

        std::cout << "Results written to " + resultDir.string() << std::endl;
}


//////////////////////////////////////////////////////////
//  RUN 3D
//////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(run3d)
{
    constexpr int dimension = 3;

    std::cout << "***********************************" << std::endl;
    std::cout << "**      Structure                **" << std::endl;
    std::cout << "***********************************" << std::endl;

    const boost::filesystem::path resultDir             (boost::filesystem::initial_path().string() + "/results_ConstraintsNodeToElement3d/");
    const boost::filesystem::path meshFilePathMatrix(boost::filesystem::initial_path().string() + "/ConstraintsNodeToElementMatrixMesh3d.msh");
    const boost::filesystem::path meshFilePathFiber(boost::filesystem::initial_path().string() + "/ConstraintsNodeToElementFiberMesh3d.msh");


    NuTo::Structure s(dimension);
    s.SetVerboseLevel(10);
    s.SetShowTime(false);

    std::cout << "***********************************" << std::endl;
    std::cout << "**      Integration Scheme       **" << std::endl;
    std::cout << "***********************************" << std::endl;

    NuTo::NewmarkDirect myIntegrationScheme(&s);
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

    auto fibreSection = NuTo::SectionTruss::Create(Parameters::mFibreCrossSection);

    std::cout << "***********************************" << std::endl;
    std::cout << "**      Material                 **" << std::endl;
    std::cout << "***********************************" << std::endl;

//    int matrixMaterial = s.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
//    s.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, Parameters::mMatrixYoungsModulus);
//    s.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, Parameters::mMatrixPoissonsRatio);

    int matrixMaterial = s.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::GRADIENT_DAMAGE_ENGINEERING_STRESS);
    s.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, Parameters::mMatrixYoungsModulus);
    s.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, Parameters::mMatrixPoissonsRatio);
    s.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::TENSILE_STRENGTH, Parameters::mMatrixTensileStrength);
    s.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::COMPRESSIVE_STRENGTH, Parameters::mMatrixCompressiveStrength);
    s.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS, Parameters::mMatrixNonlocalRadius);
    s.ConstitutiveLawSetDamageLaw(matrixMaterial, NuTo::Constitutive::DamageLawExponential::Create(Parameters::mMatrixTensileStrength / Parameters::mMatrixYoungsModulus, Parameters::mMatrixTensileStrength / Parameters::mMatrixFractureEnergy));

    int fibreMaterial = s.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    s.ConstitutiveLawSetParameterDouble(fibreMaterial, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, Parameters::mFibreYoungsModulus);
    s.ConstitutiveLawSetParameterDouble(fibreMaterial, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, Parameters::mFibrePoissonsRatio);

    std::cout << "***********************************" << std::endl;
    std::cout << "**      Interpolation Type       **" << std::endl;
    std::cout << "***********************************" << std::endl;

    int matrixInterpolationType = s.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TETRAHEDRON3D);
    s.InterpolationTypeAdd(matrixInterpolationType, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    s.InterpolationTypeAdd(matrixInterpolationType, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    s.InterpolationTypeAdd(matrixInterpolationType, NuTo::Node::eDof::NONLOCALEQSTRAIN, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);

    int fibreInterpolationType = s.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TRUSSXD);
    s.InterpolationTypeAdd(fibreInterpolationType, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    s.InterpolationTypeAdd(fibreInterpolationType, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);

    std::cout << "***********************************" << std::endl;
    std::cout << "**      Import Matrix Mesh       **" << std::endl;
    std::cout << "***********************************" << std::endl;

    auto createdGroupIdMatrix = s.ImportFromGmsh(meshFilePathMatrix.string());
    int groupIdMatrix = createdGroupIdMatrix[0].first;

    s.ElementGroupSetInterpolationType(groupIdMatrix, matrixInterpolationType);
    s.InterpolationTypeSetIntegrationType(matrixInterpolationType, NuTo::eIntegrationType::IntegrationType3D4NGauss4Ip);
    s.ElementConvertToInterpolationType(groupIdMatrix);
    s.ElementGroupSetConstitutiveLaw(groupIdMatrix, matrixMaterial);

    std::cout << "***********************************" << std::endl;
    std::cout << "**      Import Fiber Mesh        **" << std::endl;
    std::cout << "***********************************" << std::endl;

    auto createdGroupIdFiber = s.ImportFromGmsh(meshFilePathFiber.string());
    int groupIdFiber = createdGroupIdFiber[0].first;

    s.ElementGroupSetInterpolationType(groupIdFiber, fibreInterpolationType);
    s.ElementGroupSetSection(groupIdFiber, fibreSection);
    s.ElementGroupSetConstitutiveLaw(groupIdFiber, fibreMaterial);
    s.ElementConvertToInterpolationType(groupIdFiber);

    std::cout << "***********************************" << std::endl;
    std::cout << "**      Boundary Conditions      **" << std::endl;
    std::cout << "***********************************" << std::endl;

    const auto& groupNodeBCLeft = s.GroupGetNodesAtCoordinate(NuTo::eDirection::X, 0);
    s.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Constraint::Component(groupNodeBCLeft, {NuTo::eDirection::X, NuTo::eDirection::Y, NuTo::eDirection::Z}));

    const auto& nodeRight = s.NodeGetAtCoordinate(Eigen::Vector3d(10, 0, 0));
    s.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Constraint::Component(nodeRight, {NuTo::eDirection::Y, NuTo::eDirection::Z}));

    std::cout << "***********************************" << std::endl;
    std::cout << "**      Constraints              **" << std::endl;
    std::cout << "***********************************" << std::endl;


    constexpr int       numSearchDomains        = 3;
    constexpr double    length                  = 10.0;
    constexpr double    deltaLength             = length/numSearchDomains;

    for (int iDomain = 0; iDomain < numSearchDomains; ++iDomain)
    {
        int groupMatrixNodes = s.GroupCreate(NuTo::eGroupId::Nodes);
        s.GroupAddNodeFromElementGroupCoordinateRange(groupMatrixNodes, groupIdMatrix, 0, iDomain * deltaLength, (iDomain+1) * deltaLength);

        int groupMatrixElements = s.GroupCreate(NuTo::eGroupId::Elements);
        s.GroupAddElementsFromNodes(groupMatrixElements, groupMatrixNodes, false);

        int groupConstraintNodes = s.GroupCreate(NuTo::eGroupId::Nodes);
        s.GroupAddNodeFromElementGroupCoordinateRange(groupConstraintNodes, groupIdFiber, 0, iDomain * deltaLength, (iDomain+1) * deltaLength);

        for (int nodeId : s.GroupGetMemberIds(groupConstraintNodes))
            s.ConstraintLinearEquationNodeToElementCreate(nodeId, groupMatrixElements,
                                                                    NuTo::Node::eDof::DISPLACEMENTS);
    }

    std::cout << "***********************************" << std::endl;
    std::cout << "**      Loads                    **" << std::endl;
    std::cout << "***********************************" << std::endl;

    const auto& groupNodeLoadRight = s.GroupGetNodesAtCoordinate(NuTo::eDirection::X, 10);
    const auto& timeDependentConstraint = NuTo::Constraint::RhsRamp(Parameters::mSimulationTime, Parameters::mLoad);
    s.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS, 
            NuTo::Constraint::Component(groupNodeLoadRight, {NuTo::eDirection::X}, timeDependentConstraint));

    std::cout << "***********************************" << std::endl;
    std::cout << "**      Visualization            **" << std::endl;
    std::cout << "***********************************" << std::endl;

    s.AddVisualizationComponent(groupIdMatrix, NuTo::eVisualizeWhat::DISPLACEMENTS);
    s.AddVisualizationComponent(groupIdMatrix, NuTo::eVisualizeWhat::DAMAGE);
//
    s.AddVisualizationComponent(groupIdFiber, NuTo::eVisualizeWhat::DISPLACEMENTS);

    std::cout << "***********************************" << std::endl;
    std::cout << "**      Solver                   **" << std::endl;
    std::cout << "***********************************" << std::endl;

    s.CalculateMaximumIndependentSets();
    myIntegrationScheme.Solve(Parameters::mSimulationTime);

    std::cout << "***********************************" << std::endl;
    std::cout << "**      Postprocessing           **" << std::endl;
    std::cout << "***********************************" << std::endl;

    // Calculate the displacements at one of the constrained points in the matrix and the fiber and check if they match
    auto elementPtrMatrix = s.ElementGetElementPtr(0);

    // natural coordinates of (2.5, 0.20, 0.25) in tetrahedron 0
    Eigen::Vector3d nodeCoords;
    nodeCoords[0] = 0.125;
    nodeCoords[1] = 0.20;
    nodeCoords[2] = 0.125;
    auto dispInMatrix = elementPtrMatrix->InterpolateDofGlobal(nodeCoords, NuTo::Node::eDof::DISPLACEMENTS);
    auto coordsInMatrix = elementPtrMatrix->InterpolateDofGlobal(nodeCoords, NuTo::Node::eDof::COORDINATES);
    std::cout << "dispInMatrix \n" << dispInMatrix << std::endl;
    std::cout << "coordsInMatrix \n" << coordsInMatrix << std::endl;

    auto elementPtrFiber = s.ElementGetElementPtr(6);

    // natural coordinates of (2.5, 0.20, 0.25) in fiber 6
    Eigen::VectorXd fibreNaturalCoordiante = Eigen::VectorXd::Constant(1, -1);
    auto dispInFiber = elementPtrFiber->InterpolateDofGlobal(fibreNaturalCoordiante, NuTo::Node::eDof::DISPLACEMENTS);
    auto coordsInFiber = elementPtrFiber->InterpolateDofGlobal(fibreNaturalCoordiante, NuTo::Node::eDof::COORDINATES);
    std::cout << "dispInFiber \n" << dispInFiber << std::endl;
    std::cout << "coordsInFiber \n" << coordsInFiber << std::endl;

    BoostUnitTest::CheckVector(dispInFiber, dispInMatrix, 3);
    BoostUnitTest::CheckVector(coordsInFiber, coordsInMatrix, 3);

    std::cout << "Results written to " + resultDir.string() << std::endl;
}

