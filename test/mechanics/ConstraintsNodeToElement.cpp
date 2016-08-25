//============================================================================
// Name        : ConstraintsNodeToElement.cpp
// Author      : Philip Huschke
// Version     : 11 Jan 2016
// Copyright   :
// Description : Test for the constraints between nodes inside an element
//
//
//============================================================================


#include "nuto/math/FullVector.h"
#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/elements/ElementDataEnum.h"
#include "nuto/mechanics/elements/IpDataEnum.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeEnum.h"
#include "nuto/mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "nuto/mechanics/groups/GroupEnum.h"
#include "nuto/mechanics/nodes/NodeEnum.h"
#include "nuto/mechanics/sections/SectionEnum.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"
#include "nuto/visualize/VisualizeEnum.h"

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

void run2d()
{
        constexpr int dimension = 2;
        const NuTo::FullVector<double, dimension> directionX = NuTo::FullVector<double, dimension>::UnitX();
        const NuTo::FullVector<double, dimension> directionY = NuTo::FullVector<double, dimension>::UnitY();

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Structure                **" << std::endl;
        std::cout << "***********************************" << std::endl;

        const boost::filesystem::path resultDir             (boost::filesystem::initial_path().string() + "/results_ConstraintsNodeToElement2d/");
        const boost::filesystem::path meshFilePathMatrix    (boost::filesystem::initial_path().string() + "/ConstraintsNodeToElementMatrixMesh2d.msh");
        const boost::filesystem::path meshFilePathFiber     (boost::filesystem::initial_path().string() + "/ConstraintsNodeToElementFiberMesh2d.msh");

        NuTo::Structure myStructure(dimension);
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

        int matrixSection = myStructure.SectionCreate(NuTo::eSectionType::PLANE_STRESS);
        myStructure.SectionSetThickness(matrixSection, Parameters::mMatrixThickness);

        int fibreSection = myStructure.SectionCreate(NuTo::eSectionType::TRUSS);
        myStructure.SectionSetArea(fibreSection, Parameters::mFibreCrossSection);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Material                 **" << std::endl;
        std::cout << "***********************************" << std::endl;

//        int matrixMaterial = myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
//        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, Parameters::mMatrixYoungsModulus);
//        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, Parameters::mMatrixPoissonsRatio);

        int matrixMaterial = myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::GRADIENT_DAMAGE_ENGINEERING_STRESS);
        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, Parameters::mMatrixYoungsModulus);
        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, Parameters::mMatrixPoissonsRatio);
        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::TENSILE_STRENGTH, Parameters::mMatrixTensileStrength);
        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::COMPRESSIVE_STRENGTH, Parameters::mMatrixCompressiveStrength);
        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS, Parameters::mMatrixNonlocalRadius);
        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::FRACTURE_ENERGY, Parameters::mMatrixFractureEnergy);

        int fibreMaterial = myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
        myStructure.ConstitutiveLawSetParameterDouble(fibreMaterial, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, Parameters::mFibreYoungsModulus);
        myStructure.ConstitutiveLawSetParameterDouble(fibreMaterial, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, Parameters::mFibrePoissonsRatio);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Interpolation Type       **" << std::endl;
        std::cout << "***********************************" << std::endl;

        int matrixInterpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TRIANGLE2D);
        myStructure.InterpolationTypeAdd(matrixInterpolationType, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
        myStructure.InterpolationTypeAdd(matrixInterpolationType, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
        myStructure.InterpolationTypeAdd(matrixInterpolationType, NuTo::Node::eDof::NONLOCALEQSTRAIN, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);


        int fibreInterpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TRUSSXD);
        myStructure.InterpolationTypeAdd(fibreInterpolationType, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
        myStructure.InterpolationTypeAdd(fibreInterpolationType, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Import Matrix Mesh       **" << std::endl;
        std::cout << "***********************************" << std::endl;

        NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> createdGroupIdMatrix = myStructure.ImportFromGmsh(meshFilePathMatrix.string(), NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::STATICDATA);
        int groupIdMatrix = createdGroupIdMatrix.GetValue(0, 0);


        myStructure.ElementGroupSetInterpolationType(groupIdMatrix, matrixInterpolationType);
        myStructure.InterpolationTypeSetIntegrationType(matrixInterpolationType, NuTo::eIntegrationType::IntegrationType2D3NGauss3Ip, NuTo::IpData::eIpDataType::STATICDATA);
        myStructure.ElementTotalConvertToInterpolationType();
        myStructure.ElementGroupSetSection(groupIdMatrix, matrixSection);
        myStructure.ElementGroupSetConstitutiveLaw(groupIdMatrix, matrixMaterial);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Import Fiber Mesh        **" << std::endl;
        std::cout << "***********************************" << std::endl;

        NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> createdGroupIdFiber = myStructure.ImportFromGmsh(meshFilePathFiber.string(), NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::STATICDATA);
        int groupIdFiber = createdGroupIdFiber.GetValue(0, 0);

        myStructure.ElementGroupSetInterpolationType(groupIdFiber, fibreInterpolationType);
        myStructure.InterpolationTypeSetIntegrationType(fibreInterpolationType, NuTo::eIntegrationType::IntegrationType1D2NGauss3Ip, NuTo::IpData::eIpDataType::STATICDATA);
        myStructure.ElementConvertToInterpolationType(groupIdFiber);
        myStructure.ElementGroupSetSection(groupIdFiber, fibreSection);
        myStructure.ElementGroupSetConstitutiveLaw(groupIdFiber, fibreMaterial);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Boundary Conditions      **" << std::endl;
        std::cout << "***********************************" << std::endl;

        int groupNodeBCLeft = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
        myStructure.GroupAddNodeCoordinateRange(groupNodeBCLeft, 0, 0.0 - 1e-6, 0.0 + 1e-6);
        myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodeBCLeft, directionX, 0);

        NuTo::FullVector<double, dimension> nodeCoords;
        nodeCoords[0] = 0.0;
        nodeCoords[1] = 0.0;
        int nodeLeft = myStructure.NodeGetIdAtCoordinate(nodeCoords, 1e-6);
        myStructure.ConstraintLinearSetDisplacementNode(nodeLeft, directionY, 0);

        nodeCoords[0] = 10.0;
        nodeCoords[1] = 0.0;
        int nodeRight = myStructure.NodeGetIdAtCoordinate(nodeCoords, 1e-6);
        myStructure.ConstraintLinearSetDisplacementNode(nodeRight, directionY, 0);


        std::cout << "***********************************" << std::endl;
        std::cout << "**      Constraints              **" << std::endl;
        std::cout << "***********************************" << std::endl;


        constexpr int       numNearestNeighbours    = 1;
        constexpr int       numSearchDomains        = 1;
        constexpr double    length                  = 10.0;
        constexpr double    deltaLength             = length/numSearchDomains;

        for (int iDomain = 0; iDomain < numSearchDomains; ++iDomain)
        {
            int groupMatrixNodes = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
            myStructure.GroupAddNodeFromElementGroupCoordinateRange(groupMatrixNodes, groupIdMatrix, 0, iDomain * deltaLength, (iDomain+1) * deltaLength);

            int groupMatrixElements = myStructure.GroupCreate(NuTo::eGroupId::Elements);
            myStructure.GroupAddElementsFromNodes(groupMatrixElements, groupMatrixNodes, false);

            int groupConstraintNodes = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
            myStructure.GroupAddNodeFromElementGroupCoordinateRange(groupConstraintNodes, groupIdFiber, 0, iDomain * deltaLength, (iDomain+1) * deltaLength);

            auto nodeIds = myStructure.GroupGetMemberIds(groupConstraintNodes);

            for (int iNode = 0; iNode < nodeIds.rows(); ++iNode)
            {
                myStructure.ConstraintLinearEquationNodeToElementCreate(nodeIds(iNode, 0), groupMatrixElements, NuTo::Node::eDof::DISPLACEMENTS, numNearestNeighbours);
            }
        }

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Loads                    **" << std::endl;
        std::cout << "***********************************" << std::endl;

        int groupNodeLoadRight = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
        myStructure.GroupAddNodeCoordinateRange(groupNodeLoadRight, 0, 10.0 - 1e-6, 10.0 + 1e-6);

        int timeDependentConstraint = myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodeLoadRight, directionX, 1);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Visualization            **" << std::endl;
        std::cout << "***********************************" << std::endl;

        myStructure.AddVisualizationComponent(groupIdMatrix, NuTo::eVisualizeWhat::DISPLACEMENTS);
        myStructure.AddVisualizationComponent(groupIdMatrix, NuTo::eVisualizeWhat::CONSTITUTIVE);
        myStructure.AddVisualizationComponent(groupIdMatrix, NuTo::eVisualizeWhat::DAMAGE);

        myStructure.AddVisualizationComponent(groupIdFiber, NuTo::eVisualizeWhat::DISPLACEMENTS);
        myStructure.AddVisualizationComponent(groupIdFiber, NuTo::eVisualizeWhat::CONSTITUTIVE);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Solver                   **" << std::endl;
        std::cout << "***********************************" << std::endl;

        myStructure.NodeBuildGlobalDofs();
        myStructure.CalculateMaximumIndependentSets();


        NuTo::FullMatrix<double, 2, 2> timeDependentLoad;
        timeDependentLoad(0, 0) = 0;
        timeDependentLoad(1, 0) = Parameters::mSimulationTime;

        timeDependentLoad(0, 1) = 0;
        timeDependentLoad(1, 1) = Parameters::mLoad;

        myIntegrationScheme.AddTimeDependentConstraint(timeDependentConstraint, timeDependentLoad);

        myIntegrationScheme.Solve(Parameters::mSimulationTime);



        std::cout << "***********************************" << std::endl;
        std::cout << "**      Postprocessing           **" << std::endl;
        std::cout << "***********************************" << std::endl;



        // Calculate the displacements at one of the constrained points in the matrix and the fiber and check if they match
        auto elementPtrMatrix = myStructure.ElementGetElementPtr(0);

        // natural coordinates of (2.5, 0.25) in tetrahedron 0
        nodeCoords[0] = 0.5;
        nodeCoords[1] = 0.25;

        auto dispInMatrix = elementPtrMatrix->InterpolateDofGlobal(nodeCoords, NuTo::Node::eDof::DISPLACEMENTS);
        auto coordsInMatrix = elementPtrMatrix->InterpolateDofGlobal(nodeCoords, NuTo::Node::eDof::COORDINATES);
        std::cout << "dispInMatrix \n" << dispInMatrix << std::endl;
        std::cout << "coordsInMatrix \n" << coordsInMatrix << std::endl;

        auto elementPtrFiber = myStructure.ElementGetElementPtr(4);

        // natural coordinates of (2.5, 0.25) in fiber 4
        nodeCoords[0] = -1.0;
        nodeCoords[1] = 0.0;
        auto dispInFiber = elementPtrFiber->InterpolateDofGlobal(nodeCoords, NuTo::Node::eDof::DISPLACEMENTS);
        auto coordsInFiber = elementPtrFiber->InterpolateDofGlobal(nodeCoords, NuTo::Node::eDof::COORDINATES);
        std::cout << "dispInFiber \n" << dispInFiber << std::endl;
        std::cout << "coordsInFiber \n" << coordsInFiber << std::endl;


        if ( (dispInMatrix - dispInFiber).norm() > 1e-6 or (coordsInMatrix - coordsInFiber).norm() > 1e-6)
            throw NuTo::MechanicsException(std::string(__PRETTY_FUNCTION__) + ": \t Displacements and/or coordinates of fiber and matrix do not match!");


        std::cout << "Results written to " + resultDir.string() << std::endl;

}


//////////////////////////////////////////////////////////
//  RUN 3D
//////////////////////////////////////////////////////////

void run3d()
{
    constexpr int dimension = 3;
    const NuTo::FullVector<double, dimension> directionX = NuTo::FullVector<double, dimension>::UnitX();
    const NuTo::FullVector<double, dimension> directionY = NuTo::FullVector<double, dimension>::UnitY();
    const NuTo::FullVector<double, dimension> directionZ = NuTo::FullVector<double, dimension>::UnitZ();

    std::cout << "***********************************" << std::endl;
    std::cout << "**      Structure                **" << std::endl;
    std::cout << "***********************************" << std::endl;

    const boost::filesystem::path resultDir             (boost::filesystem::initial_path().string() + "/results_ConstraintsNodeToElement3d/");
    const boost::filesystem::path meshFilePathMatrix(boost::filesystem::initial_path().string() + "/ConstraintsNodeToElementMatrixMesh3d.msh");
    const boost::filesystem::path meshFilePathFiber(boost::filesystem::initial_path().string() + "/ConstraintsNodeToElementFiberMesh3d.msh");


    NuTo::Structure myStructure(dimension);
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

    int matrixSection = myStructure.SectionCreate(NuTo::eSectionType::VOLUME);

    int fibreSection = myStructure.SectionCreate(NuTo::eSectionType::TRUSS);
    myStructure.SectionSetArea(fibreSection, Parameters::mFibreCrossSection);

    std::cout << "***********************************" << std::endl;
    std::cout << "**      Material                 **" << std::endl;
    std::cout << "***********************************" << std::endl;

//    int matrixMaterial = myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
//    myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, Parameters::mMatrixYoungsModulus);
//    myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, Parameters::mMatrixPoissonsRatio);

    int matrixMaterial = myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::GRADIENT_DAMAGE_ENGINEERING_STRESS);
    myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, Parameters::mMatrixYoungsModulus);
    myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, Parameters::mMatrixPoissonsRatio);
    myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::TENSILE_STRENGTH, Parameters::mMatrixTensileStrength);
    myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::COMPRESSIVE_STRENGTH, Parameters::mMatrixCompressiveStrength);
    myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS, Parameters::mMatrixNonlocalRadius);
    myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::FRACTURE_ENERGY, Parameters::mMatrixFractureEnergy);

    int fibreMaterial = myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    myStructure.ConstitutiveLawSetParameterDouble(fibreMaterial, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, Parameters::mFibreYoungsModulus);
    myStructure.ConstitutiveLawSetParameterDouble(fibreMaterial, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, Parameters::mFibrePoissonsRatio);

    std::cout << "***********************************" << std::endl;
    std::cout << "**      Interpolation Type       **" << std::endl;
    std::cout << "***********************************" << std::endl;

    int matrixInterpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TETRAHEDRON3D);
    myStructure.InterpolationTypeAdd(matrixInterpolationType, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(matrixInterpolationType, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    myStructure.InterpolationTypeAdd(matrixInterpolationType, NuTo::Node::eDof::NONLOCALEQSTRAIN, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);

    int fibreInterpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TRUSSXD);
    myStructure.InterpolationTypeAdd(fibreInterpolationType, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(fibreInterpolationType, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);

    std::cout << "***********************************" << std::endl;
    std::cout << "**      Import Matrix Mesh       **" << std::endl;
    std::cout << "***********************************" << std::endl;

    NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> createdGroupIdMatrix = myStructure.ImportFromGmsh(meshFilePathMatrix.string(), NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::STATICDATA);
    int groupIdMatrix = createdGroupIdMatrix.GetValue(0, 0);

    myStructure.ElementGroupSetInterpolationType(groupIdMatrix, matrixInterpolationType);
    myStructure.InterpolationTypeSetIntegrationType(matrixInterpolationType, NuTo::eIntegrationType::IntegrationType3D4NGauss4Ip, NuTo::IpData::eIpDataType::STATICDATA);
    myStructure.ElementConvertToInterpolationType(groupIdMatrix);
    myStructure.ElementGroupSetSection(groupIdMatrix, matrixSection);
    myStructure.ElementGroupSetConstitutiveLaw(groupIdMatrix, matrixMaterial);

    std::cout << "***********************************" << std::endl;
    std::cout << "**      Import Fiber Mesh        **" << std::endl;
    std::cout << "***********************************" << std::endl;

    NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> createdGroupIdFiber = myStructure.ImportFromGmsh(meshFilePathFiber.string(), NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::STATICDATA);
    int groupIdFiber = createdGroupIdFiber.GetValue(0, 0);

    myStructure.ElementGroupSetInterpolationType(groupIdFiber, fibreInterpolationType);
    myStructure.ElementGroupSetSection(groupIdFiber, fibreSection);
    myStructure.ElementGroupSetConstitutiveLaw(groupIdFiber, fibreMaterial);
    myStructure.ElementConvertToInterpolationType(groupIdFiber);

    std::cout << "***********************************" << std::endl;
    std::cout << "**      Boundary Conditions      **" << std::endl;
    std::cout << "***********************************" << std::endl;

    int groupNodeBCLeft = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    myStructure.GroupAddNodeCoordinateRange(groupNodeBCLeft, 0, 0.0 - 1e-6, 0.0 + 1e-6);
    myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodeBCLeft, directionX, 0);
    myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodeBCLeft, directionY, 0);
    myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodeBCLeft, directionZ, 0);


    NuTo::FullVector<double, dimension> nodeCoords;

    nodeCoords[0] = 10.0;
    nodeCoords[1] = 0.0;
    nodeCoords[2] = 0.0;
    int nodeRight = myStructure.NodeGetIdAtCoordinate(nodeCoords, 1e-6);
    myStructure.ConstraintLinearSetDisplacementNode(nodeRight, directionY, 0);
    myStructure.ConstraintLinearSetDisplacementNode(nodeRight, directionZ, 0);

    std::cout << "***********************************" << std::endl;
    std::cout << "**      Constraints              **" << std::endl;
    std::cout << "***********************************" << std::endl;


    constexpr int       numNearestNeighbours    = 1;
    constexpr int       numSearchDomains        = 3;
    constexpr double    length                  = 10.0;
    constexpr double    deltaLength             = length/numSearchDomains;

    for (int iDomain = 0; iDomain < numSearchDomains; ++iDomain)
    {
        int groupMatrixNodes = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
        myStructure.GroupAddNodeFromElementGroupCoordinateRange(groupMatrixNodes, groupIdMatrix, 0, iDomain * deltaLength, (iDomain+1) * deltaLength);

        int groupMatrixElements = myStructure.GroupCreate(NuTo::eGroupId::Elements);
        myStructure.GroupAddElementsFromNodes(groupMatrixElements, groupMatrixNodes, false);

        int groupConstraintNodes = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
        myStructure.GroupAddNodeFromElementGroupCoordinateRange(groupConstraintNodes, groupIdFiber, 0, iDomain * deltaLength, (iDomain+1) * deltaLength);

        auto nodeIds = myStructure.GroupGetMemberIds(groupConstraintNodes);

        for (int iNode = 0; iNode < nodeIds.rows(); ++iNode)
        {
            myStructure.ConstraintLinearEquationNodeToElementCreate(nodeIds(iNode, 0), groupMatrixElements, NuTo::Node::eDof::DISPLACEMENTS, numNearestNeighbours);
        }
    }

    std::cout << "***********************************" << std::endl;
    std::cout << "**      Loads                    **" << std::endl;
    std::cout << "***********************************" << std::endl;

    int groupNodeLoadRight = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    myStructure.GroupAddNodeCoordinateRange(groupNodeLoadRight, 0, 10.0 - 1e-6, 10.0 + 1e-6);

    int timeDependentConstraint = myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodeLoadRight, directionX, 1);

    std::cout << "***********************************" << std::endl;
    std::cout << "**      Visualization            **" << std::endl;
    std::cout << "***********************************" << std::endl;

    myStructure.AddVisualizationComponent(groupIdMatrix, NuTo::eVisualizeWhat::DISPLACEMENTS);
    myStructure.AddVisualizationComponent(groupIdMatrix, NuTo::eVisualizeWhat::CONSTITUTIVE);
    myStructure.AddVisualizationComponent(groupIdMatrix, NuTo::eVisualizeWhat::ELEMENT);
    myStructure.AddVisualizationComponent(groupIdMatrix, NuTo::eVisualizeWhat::DAMAGE);
//
    myStructure.AddVisualizationComponent(groupIdFiber, NuTo::eVisualizeWhat::DISPLACEMENTS);
    myStructure.AddVisualizationComponent(groupIdFiber, NuTo::eVisualizeWhat::CONSTITUTIVE);
    myStructure.AddVisualizationComponent(groupIdFiber, NuTo::eVisualizeWhat::ELEMENT);

    std::cout << "***********************************" << std::endl;
    std::cout << "**      Solver                   **" << std::endl;
    std::cout << "***********************************" << std::endl;

    myStructure.NodeBuildGlobalDofs();
    myStructure.CalculateMaximumIndependentSets();


    NuTo::FullMatrix<double, 2, 2> timeDependentLoad;
    timeDependentLoad(0, 0) = 0;
    timeDependentLoad(1, 0) = Parameters::mSimulationTime;

    timeDependentLoad(0, 1) = 0;
    timeDependentLoad(1, 1) = Parameters::mLoad;

    myIntegrationScheme.AddTimeDependentConstraint(timeDependentConstraint, timeDependentLoad);

    myIntegrationScheme.Solve(Parameters::mSimulationTime);



    std::cout << "***********************************" << std::endl;
    std::cout << "**      Postprocessing           **" << std::endl;
    std::cout << "***********************************" << std::endl;

    // Calculate the displacements at one of the constrained points in the matrix and the fiber and check if they match
    auto elementPtrMatrix = myStructure.ElementGetElementPtr(0);

    // natural coordinates of (2.5, 0.20, 0.25) in tetrahedron 0
    nodeCoords[0] = 0.125;
    nodeCoords[1] = 0.20;
    nodeCoords[2] = 0.125;
    auto dispInMatrix = elementPtrMatrix->InterpolateDofGlobal(nodeCoords, NuTo::Node::eDof::DISPLACEMENTS);
    auto coordsInMatrix = elementPtrMatrix->InterpolateDofGlobal(nodeCoords, NuTo::Node::eDof::COORDINATES);
    std::cout << "dispInMatrix \n" << dispInMatrix << std::endl;
    std::cout << "coordsInMatrix \n" << coordsInMatrix << std::endl;

    auto elementPtrFiber = myStructure.ElementGetElementPtr(6);

    // natural coordinates of (2.5, 0.20, 0.25) in fiber 6
    nodeCoords[0] = -1.0;
    nodeCoords[1] = 0.0;
    nodeCoords[2] = 0.0;
    auto dispInFiber = elementPtrFiber->InterpolateDofGlobal(nodeCoords, NuTo::Node::eDof::DISPLACEMENTS);
    auto coordsInFiber = elementPtrFiber->InterpolateDofGlobal(nodeCoords, NuTo::Node::eDof::COORDINATES);
    std::cout << "dispInFiber \n" << dispInFiber << std::endl;
    std::cout << "coordsInFiber \n" << coordsInFiber << std::endl;

    if ( (dispInMatrix - dispInFiber).norm() > 1e-6 or (coordsInMatrix - coordsInFiber).norm() > 1e-6)
        throw NuTo::MechanicsException(std::string(__PRETTY_FUNCTION__) + ": \t Displacements and/or coordinates of fiber and matrix do not match!");


    std::cout << "Results written to " + resultDir.string() << std::endl;


}


//////////////////////////////////////////////////////////
//  MAIN
//////////////////////////////////////////////////////////

int main()
{
    try
    {
        run2d();
        run3d();

    } catch (NuTo::MechanicsException& e)
    {
        std::cout << e.ErrorMessage() << std::endl;
        return EXIT_FAILURE;

    } catch (NuTo::MathException& e)
    {
        std::cout << e.ErrorMessage() << std::endl;
        return EXIT_FAILURE;

    } catch (...)
    {
        std::cout << "Something else went wrong :-(" << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "***********************************" << std::endl;
    std::cout << "**      End                      **" << std::endl;
    std::cout << "***********************************" << std::endl;

}

