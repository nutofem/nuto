/*
 * Interpolation2D.cpp
 *
 *  Created on: 20 Mar 2015
 *      Author: ttitsche
 */

#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/elements/ElementBase.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/structures/StructureOutputBlockMatrix.h"
#include "mechanics/structures/StructureOutputBlockVector.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/interpolationtypes/InterpolationType.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "mechanics/interpolationtypes/Interpolation2DTriangle.h"


#include "mechanics/integrationtypes/IntegrationType1D2NGauss2Ip.h"
#include "mechanics/integrationtypes/IntegrationType1D2NGauss3Ip.h"
#include "mechanics/integrationtypes/IntegrationType1D2NGauss4Ip.h"

#include "mechanics/integrationtypes/IntegrationType2D3NGauss13Ip.h"
#include "mechanics/integrationtypes/IntegrationType2D4NGauss4Ip.h"
#include "mechanics/integrationtypes/IntegrationType3D4NGauss4Ip.h"
#include "mechanics/integrationtypes/IntegrationType3D8NGauss2x2x2Ip.h"
#include "mechanics/integrationtypes/IntegrationType3D6NGauss2x3Ip.h"
#include "mechanics/integrationtypes/IntegrationTypeEnum.h"

#include "base/Exception.h"
#include <boost/filesystem.hpp>

#include <map>
#include <string>
#include <iostream>

#include <cmath>

const bool PRINTRESULT = true;

//! @brief checks if the shape functions sum up to 1 (at all integration points)
void CheckPartitionOfUnity(const NuTo::InterpolationType &rIT, const NuTo::Node::eDof &dofType)
{
    for (int iIP = 0; iIP < rIT.GetCurrentIntegrationType().GetNumIntegrationPoints(); ++iIP)
    {
        Eigen::VectorXd ip = rIT.GetCurrentIntegrationType().GetLocalIntegrationPointCoordinates(iIP);
        if (abs(rIT.Get(dofType).CalculateShapeFunctions(ip).sum() - 1.) > 1.e-6)
            throw NuTo::MechanicsException("[CheckPartitionOfUnity] shape functions at integration point "
                                               + std::to_string(iIP) + " do not sum up to 1.");
    }
}

//! @brief checks wheather B = dN/dxi around node 0
void CheckDerivatives(NuTo::InterpolationType& rIT)
{
    auto& IT = rIT.Get(NuTo::Node::eDof::COORDINATES);

    for (int i = 0; i < rIT.GetNumNodes(); ++i)
    {
        auto nodeCoordinates = rIT.GetNaturalNodeCoordinates(i);

        auto B = IT.CalculateDerivativeShapeFunctionsNatural(nodeCoordinates);
        auto N = IT.CalculateShapeFunctions(nodeCoordinates);

        auto B_CDF = B;
        B_CDF.setZero();

        double delta = 1.e-6;
        for (int iDim = 0; iDim < nodeCoordinates.rows(); ++iDim)
        {
            nodeCoordinates[iDim] += delta;
            B_CDF.col(iDim) = (IT.CalculateShapeFunctions(nodeCoordinates) - N) / delta;
            nodeCoordinates[iDim] -= delta;
        }

        if ((B - B_CDF).cwiseAbs().maxCoeff() > 1.e-4)
        {
            std::cout << "B\n" << B << std::endl;
            std::cout << "B_CDF\n" << B_CDF << std::endl;
            throw NuTo::MechanicsException("[CheckDerivatives] B != dN/dXi");
        }
    }
    std::cout << "[CheckDerivatives] OK!" << std::endl;
}

//! @brief checks, whether or not the natural node coordinates match the shape functions
//! This should be true: N_i(xi_j) == 1 for i == j    and      N_j(xi_i) == 0 for i != j
void CheckShapeFunctionsAndNodePositions(NuTo::InterpolationType& rIT, int rNumNodesExpected)
{
    auto dofType = NuTo::Node::eDof::COORDINATES;
    assert(rIT.IsDof(dofType)); // coordinates exist?

    if (rIT.GetNumNodes() != rNumNodesExpected)
    {
        std::cout << rIT.GetNumNodes() << std::endl;
        std::cout << rIT.Info() << std::endl;
        throw NuTo::MechanicsException("[CheckShapeFunctionsAndNodePositions] Wrong node number");
    }


    int numNodes = rIT.GetNumNodes();
    for (int iNode = 0; iNode < numNodes; ++iNode)
    {
        auto naturalNodeCoordinate = rIT.Get(dofType).CalculateNaturalNodeCoordinates(iNode);
        auto shapeFunctions = rIT.Get(dofType).CalculateShapeFunctions(naturalNodeCoordinate);
        int numShapeFunctions = shapeFunctions.rows();
        if (numShapeFunctions != numNodes)
            throw NuTo::MechanicsException("[CheckShapeFunctionsAndNodePositions] number of nodes and number of shape functions does not match.");

        for (int iShapeFunctions = 0; iShapeFunctions < numShapeFunctions; ++iShapeFunctions)
        {
            if (iShapeFunctions == iNode)
            {
                // should be 1
                if (std::abs(shapeFunctions(iShapeFunctions) - 1) > 1.e-10)
                {
                    std::cout << "Error at ShapeFunction " << iShapeFunctions << std::endl;
                    throw NuTo::MechanicsException(
                        "[CheckShapeFunctionsAndNodePositions] shape functions and node positions do not match (should be 1).");
                }
            }
            else
            {
                // should be 0
                if (std::abs(shapeFunctions(iShapeFunctions)) > 1.e-10)
                {
                    std::cout << "Error at ShapeFunction " << iShapeFunctions << " and node " << iNode << std::endl;
                    throw NuTo::MechanicsException(
                        "[CheckShapeFunctionsAndNodePositions] shape functions and node positions do not match (should be 0).");
                }
            }
        }
    }
    std::cout << "[CheckShapeFunctionsAndNodePositions] OK!" << std::endl;

    CheckPartitionOfUnity(rIT, dofType);
    CheckDerivatives(rIT);
}


//! @brief the global index (0..mNumNodesTotal) of the node is compared to the
//! dof index (0..mNumNodes(dof)) via the local node positions
void CheckNodeIndexing(NuTo::InterpolationType& rIT)
{
    for (auto dofType : rIT.GetDofs())
    {

        int numNodesDof = rIT.Get(dofType).GetNumNodes();
        for (int iNodeDof = 0; iNodeDof < numNodesDof; ++iNodeDof)
        {
            int globalNodeIndex = rIT.Get(dofType).GetNodeIndex(iNodeDof);

            auto coordinatesLocal = rIT.Get(dofType).CalculateNaturalNodeCoordinates(iNodeDof);
            auto coordinatesGlobal = rIT.GetNaturalNodeCoordinates(globalNodeIndex);

            bool coordinatesAreEqual = true;
            for (int iDim = 0; iDim < coordinatesGlobal.rows(); ++iDim)
                if (std::abs(coordinatesGlobal(iDim) - coordinatesLocal(iDim)) > 1.e-10)
                    coordinatesAreEqual = false;

            if (not coordinatesAreEqual)
                throw NuTo::MechanicsException("[CheckNodeIndexing] Node indexing wrong.");
        }
    }
    std::cout << "[CheckNodeIndexing] OK!" << std::endl;
}

void CheckTruss()
{
        NuTo::IntegrationType1D2NGauss2Ip myIntegrationType;
    {
        NuTo::InterpolationType myIT(NuTo::Interpolation::eShapeType::TRUSS1D, 1);
        myIT.AddDofInterpolation(NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
        myIT.UpdateIntegrationType(myIntegrationType);
        CheckShapeFunctionsAndNodePositions(myIT, 3);
        myIT.PrintNodeCoordinates();
        myIT.PrintNodeIndices();
        CheckNodeIndexing (myIT);
    }


    {
        NuTo::InterpolationType myIT(NuTo::Interpolation::eShapeType::TRUSS1D, 1);
        myIT.AddDofInterpolation(NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT3);
        myIT.UpdateIntegrationType(myIntegrationType);
        CheckShapeFunctionsAndNodePositions(myIT, 4);
        myIT.PrintNodeCoordinates();
        myIT.PrintNodeIndices();
        CheckNodeIndexing (myIT);
    }


    {
        NuTo::InterpolationType myIT(NuTo::Interpolation::eShapeType::TRUSS1D, 1);
        myIT.AddDofInterpolation(NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
        myIT.AddDofInterpolation(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT4);
        myIT.AddDofInterpolation(NuTo::Node::eDof::NONLOCALEQSTRAIN, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
        myIT.UpdateIntegrationType(myIntegrationType);
        myIT.PrintNodeCoordinates();
        myIT.PrintNodeIndices();
        CheckNodeIndexing (myIT);
    }


    std::cout << "[CheckTruss] OK!" << std::endl;

}

void CheckTriangle()
{
    NuTo::IntegrationType2D3NGauss13Ip myIntegrationType;

    NuTo::InterpolationType myIT3(NuTo::Interpolation::eShapeType::TRIANGLE2D, 2);
    myIT3.AddDofInterpolation(NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    myIT3.UpdateIntegrationType(myIntegrationType);
    CheckShapeFunctionsAndNodePositions(myIT3, 3);

    NuTo::InterpolationType myIT6(NuTo::Interpolation::eShapeType::TRIANGLE2D, 2);
    myIT6.AddDofInterpolation(NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    myIT6.UpdateIntegrationType(myIntegrationType);
    CheckShapeFunctionsAndNodePositions(myIT6, 6);

    NuTo::InterpolationType myIT10(NuTo::Interpolation::eShapeType::TRIANGLE2D, 2);
    myIT10.AddDofInterpolation(NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT3);
    myIT10.UpdateIntegrationType(myIntegrationType);
    CheckShapeFunctionsAndNodePositions(myIT10, 10);

    NuTo::InterpolationType myIT15(NuTo::Interpolation::eShapeType::TRIANGLE2D, 2);
    myIT15.AddDofInterpolation(NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT4);
    myIT15.UpdateIntegrationType(myIntegrationType);
    CheckShapeFunctionsAndNodePositions(myIT15, 15);

    NuTo::InterpolationType myIT(NuTo::Interpolation::eShapeType::TRIANGLE2D, 2);
    myIT.AddDofInterpolation(NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    myIT.AddDofInterpolation(NuTo::Node::eDof::TEMPERATURE, NuTo::Interpolation::eTypeOrder::EQUIDISTANT4);
    myIT.AddDofInterpolation(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT3);
    myIT.AddDofInterpolation(NuTo::Node::eDof::NONLOCALEQSTRAIN, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    myIT.UpdateIntegrationType(myIntegrationType);

    if (PRINTRESULT)
    {
        myIT.PrintNodeCoordinates();
        myIT.PrintNodeIndices();
    }

    CheckNodeIndexing(myIT);
}

void CheckQuad()
{

    NuTo::IntegrationType2D4NGauss4Ip myIntegrationType;

    NuTo::InterpolationType myIT4(NuTo::Interpolation::eShapeType::QUAD2D, 2);
    myIT4.AddDofInterpolation(NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    myIT4.UpdateIntegrationType(myIntegrationType);
    CheckShapeFunctionsAndNodePositions(myIT4, 4);

    NuTo::InterpolationType myIT8(NuTo::Interpolation::eShapeType::QUAD2D, 2);
    myIT8.AddDofInterpolation(NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    myIT8.UpdateIntegrationType(myIntegrationType);
    CheckShapeFunctionsAndNodePositions(myIT8, 8);

    NuTo::InterpolationType myIT9(NuTo::Interpolation::eShapeType::QUAD2D, 2);
    myIT9.AddDofInterpolation(NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::LOBATTO2);
    myIT9.UpdateIntegrationType(myIntegrationType);
    CheckShapeFunctionsAndNodePositions(myIT9, 9);

    NuTo::InterpolationType myIT16(NuTo::Interpolation::eShapeType::QUAD2D, 2);
    myIT16.AddDofInterpolation(NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::LOBATTO3);
    myIT16.UpdateIntegrationType(myIntegrationType);
    CheckShapeFunctionsAndNodePositions(myIT16, 16);

    NuTo::InterpolationType myIT25(NuTo::Interpolation::eShapeType::QUAD2D, 2);
    myIT25.AddDofInterpolation(NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::LOBATTO4);
    myIT25.UpdateIntegrationType(myIntegrationType);
    CheckShapeFunctionsAndNodePositions(myIT25, 25);

}

void CheckTetrahedron()
{
    NuTo::IntegrationType3D4NGauss4Ip myIntegrationType;

    NuTo::InterpolationType myIT4(NuTo::Interpolation::eShapeType::TETRAHEDRON3D, 3);
    myIT4.AddDofInterpolation(NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    myIT4.UpdateIntegrationType(myIntegrationType);
    CheckShapeFunctionsAndNodePositions(myIT4, 4);

    NuTo::InterpolationType myIT10(NuTo::Interpolation::eShapeType::TETRAHEDRON3D, 3);
    myIT10.AddDofInterpolation(NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    myIT10.UpdateIntegrationType(myIntegrationType);
    CheckShapeFunctionsAndNodePositions(myIT10, 10);

}

void CheckBrick()
{
    NuTo::IntegrationType3D8NGauss2x2x2Ip myIntegrationType;

    NuTo::InterpolationType myIT8(NuTo::Interpolation::eShapeType::BRICK3D, 3);
    myIT8.AddDofInterpolation(NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    myIT8.UpdateIntegrationType(myIntegrationType);
    CheckShapeFunctionsAndNodePositions(myIT8, 8);

    NuTo::InterpolationType myIT20(NuTo::Interpolation::eShapeType::BRICK3D, 3);
    myIT20.AddDofInterpolation(NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    myIT20.UpdateIntegrationType(myIntegrationType);
    CheckShapeFunctionsAndNodePositions(myIT20, 20);

    NuTo::InterpolationType myIT27(NuTo::Interpolation::eShapeType::BRICK3D, 3);
    myIT27.AddDofInterpolation(NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::LOBATTO2);
    myIT27.UpdateIntegrationType(myIntegrationType);
    CheckShapeFunctionsAndNodePositions(myIT27, 27);

    NuTo::InterpolationType myIT64(NuTo::Interpolation::eShapeType::BRICK3D, 3);
    myIT64.AddDofInterpolation(NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::LOBATTO3);
    myIT64.UpdateIntegrationType(myIntegrationType);
    CheckShapeFunctionsAndNodePositions(myIT64, 64);

    NuTo::InterpolationType myIT125(NuTo::Interpolation::eShapeType::BRICK3D, 3);
    myIT125.AddDofInterpolation(NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::LOBATTO4);
    myIT125.UpdateIntegrationType(myIntegrationType);
    CheckShapeFunctionsAndNodePositions(myIT125, 125);
}

void CheckPrism()
{
    NuTo::IntegrationType3D6NGauss2x3Ip myIntegrationType;

    NuTo::InterpolationType myIT6(NuTo::Interpolation::eShapeType::PRISM3D, 3);
    myIT6.AddDofInterpolation(NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    myIT6.UpdateIntegrationType(myIntegrationType);
    CheckShapeFunctionsAndNodePositions(myIT6, 6);

    NuTo::InterpolationType myIT15(NuTo::Interpolation::eShapeType::PRISM3D, 3);
    myIT15.AddDofInterpolation(NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    myIT15.UpdateIntegrationType(myIntegrationType);
    CheckShapeFunctionsAndNodePositions(myIT15, 18);
}

//! @brief API of the interpolation types
void CheckAPI()
{
    NuTo::Structure myStructure(2);

    double lx = 3;
    double ly = 5;

    //create constitutive law
    int myMatLin = myStructure.ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress");
    myStructure.ConstitutiveLawSetParameterDouble(myMatLin, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, 10);
    myStructure.ConstitutiveLawSetParameterDouble(myMatLin, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, 0.2);

    //create section
    int mySection = myStructure.SectionCreate("Plane_Strain");
    myStructure.SectionSetThickness(mySection, 1);


    // nodes without dofs, only coordinates
    int nodeIndex0 = myStructure.NodeCreate(Eigen::Vector2d({0,0}));
    int nodeIndex1 = myStructure.NodeCreate(Eigen::Vector2d({lx,0}));
    int nodeIndex2 = myStructure.NodeCreate(Eigen::Vector2d({0,ly}));

    std::vector<int> nodeIndices({nodeIndex0, nodeIndex1, nodeIndex2});

    // create interpolation type
    int myInterpolationTypeIndex = myStructure.InterpolationTypeCreate("TRIANGLE2D");

    myStructure.InterpolationTypeAdd(myInterpolationTypeIndex, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(myInterpolationTypeIndex, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    myStructure.InterpolationTypeSetIntegrationType(myInterpolationTypeIndex, NuTo::eIntegrationType::IntegrationType2D3NGauss13Ip);

    int elementIndex = myStructure.ElementCreate(myInterpolationTypeIndex, nodeIndices);
    myStructure.ElementSetConstitutiveLaw(elementIndex, myMatLin);

    int elementGroup = myStructure.GroupCreate("ELEMENTS");
    myStructure.GroupAddElement(elementGroup, elementIndex);

    myStructure.ElementConvertToInterpolationType(elementGroup, 0.01, 2.);

    int numIP = myStructure.ElementGetElementPtr(elementIndex)->GetNumIntegrationPoints();
    std::cout << numIP << std::endl;
    assert(numIP == 13);

    //assign constitutive law
    myStructure.ElementTotalSetSection(mySection);

    myStructure.NodeBuildGlobalDofs();
    myStructure.CalculateMaximumIndependentSets();

    auto hessian = myStructure.BuildGlobalHessian0();
    auto internalGradient = myStructure.BuildGlobalInternalGradient();

//    hessian.JJ.ExportToFullMatrix().Info(10,3, true);
//    internalGradient.J.Export().Info(10,3, true);

}

//! @brief Imports a mesh file and builds the hessian and the internal gradient.
void ImportFromGmsh(std::string rMeshFile)
{
    NuTo::Structure myStructure(2);
    auto groupIndices = myStructure.ImportFromGmsh(rMeshFile);

    std::cout << groupIndices.size() << std::endl;

    int interpolationType = groupIndices[0].second;
    myStructure.InterpolationTypeAdd(interpolationType, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);

    myStructure.SetVerboseLevel(10);
    myStructure.ElementConvertToInterpolationType(groupIndices[0].first);

    myStructure.InterpolationTypeSetIntegrationType(interpolationType, NuTo::eIntegrationType::IntegrationType2D3NGauss3Ip);

    myStructure.InterpolationTypeInfo(0);

    myStructure.NodeBuildGlobalDofs();
    int section = myStructure.SectionCreate("PLANE_STRAIN");
    myStructure.SectionSetThickness(section, 3);

    myStructure.ElementTotalSetSection(section);

    int constitutiveLaw = myStructure.ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress");
    myStructure.ConstitutiveLawSetParameterDouble(constitutiveLaw, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, 42);
    myStructure.ConstitutiveLawSetParameterDouble(constitutiveLaw, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, .25);
    myStructure.ElementTotalSetConstitutiveLaw(constitutiveLaw);

    myStructure.CalculateMaximumIndependentSets();
    auto hessian = myStructure.BuildGlobalHessian0();
    auto internalGradient = myStructure.BuildGlobalInternalGradient();
}


//! @brief Creates two elements from the same set of nodes. One with positive jacobian, one with a negative one. In the latter one, the nodes should be swapped.
void NodeReordering()
{
    bool errorsOccurred = false;

    // structures
    NuTo::Structure myStructureTruss(1);
    NuTo::Structure myStructureTriangle(2);
    NuTo::Structure myStructureQuad(2);
    NuTo::Structure myStructureTetrahedron(3);
    NuTo::Structure myStructureBrick(3);

    // interpolation types
    int itTruss       = myStructureTruss.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TRUSS1D);
    int itTriangle    = myStructureTriangle.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TRIANGLE2D);
    int itQuad        = myStructureQuad.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::QUAD2D);
    int itTetrahedron = myStructureTetrahedron.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TETRAHEDRON3D);
    int itBricks      = myStructureBrick.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::BRICK3D);

    myStructureTruss.InterpolationTypeAdd(itTruss, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    myStructureTriangle.InterpolationTypeAdd(itTriangle, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    myStructureQuad.InterpolationTypeAdd(itQuad, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    myStructureTetrahedron.InterpolationTypeAdd(itTetrahedron, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    myStructureBrick.InterpolationTypeAdd(itBricks, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);




    // **************
    // **  TRUSS
    // **************
    try
    {
        Eigen::MatrixXd nodes(1,2);
        nodes << 0, 1;

        std::vector<int> nodeIds = myStructureTruss.NodesCreate(nodes);
        std::vector<int> ids(2);

        // right numbering
        ids = {nodeIds[0], nodeIds[1]};
        myStructureTruss.ElementCreate(itTruss, ids);

        // "wrong" numbering
        ids = {nodeIds[1], nodeIds[0]};
        myStructureTruss.ElementCreate(itTruss, ids);
    }
    catch (NuTo::MechanicsException& e)
    {
        std::cout << "Truss failed." << std::endl;
        std::cout << e.ErrorMessage() << std::endl;
        errorsOccurred = true;
    }

    // **************
    // **  TRIANGLE
    // **************
    try
    {

        Eigen::MatrixXd nodes(2,3);
        nodes << 0, 2, 2,
                 0, 0, 3;

        std::vector<int> nodeIds = myStructureTriangle.NodesCreate(nodes);
        std::vector<int> ids(3);

        // right numbering
        ids = {nodeIds[0], nodeIds[1], nodeIds[2]};
        myStructureTriangle.ElementCreate(itTriangle, ids);

        // "wrong" numbering
        ids = {nodeIds[2], nodeIds[1], nodeIds[0]};
        myStructureTriangle.ElementCreate(itTriangle, ids);

    }
    catch (NuTo::MechanicsException& e)
    {
        std::cout << "Triangle failed." << std::endl;
        std::cout << e.ErrorMessage() << std::endl;
        errorsOccurred = true;
    }

    // **************
    // **    QUAD
    // **************
    try
    {
        Eigen::MatrixXd nodes(2,4);
        nodes << 0, 2, 2, 0,
                 0, 0, 3, 2;

        std::vector<int> nodeIds = myStructureQuad.NodesCreate(nodes);
        std::vector<int> ids(4);

        // right numbering
        ids = {nodeIds[0], nodeIds[1], nodeIds[2], nodeIds[3]};
        myStructureQuad.ElementCreate(itQuad, ids);

        // "wrong" numbering
        ids = {nodeIds[3], nodeIds[2], nodeIds[1], nodeIds[0]};
        myStructureQuad.ElementCreate(itQuad, ids);
    }
    catch (NuTo::MechanicsException& e)
    {
        std::cout << "Quad failed." << std::endl;
        std::cout << e.ErrorMessage() << std::endl;
        errorsOccurred = true;
    }

    // **************
    // ** TETRAHEDRON
    // **************
    try
    {
        Eigen::MatrixXd nodes(3,4);
        nodes << 0, 2, 0, 0,
                 0, 0, 3, 0,
                 0, 0, 0, 4;


        std::vector<int> nodeIds = myStructureTetrahedron.NodesCreate(nodes);
        std::vector<int> ids(4);

        // right numbering
        ids = {nodeIds[0], nodeIds[1], nodeIds[2], nodeIds[3]};
        myStructureTetrahedron.ElementCreate(itTetrahedron, ids);

        // "wrong" numbering
        ids = {nodeIds[0], nodeIds[3], nodeIds[2], nodeIds[1]};
        myStructureTetrahedron.ElementCreate(itTetrahedron, ids);
    }
    catch (NuTo::MechanicsException& e)
    {
        std::cout << "Tetrahedron failed." << std::endl;
        std::cout << e.ErrorMessage() << std::endl;
        errorsOccurred = true;
    }

    // **************
    // **  BRICK
    // **************
    try
    {
        Eigen::MatrixXd nodes(3,8);
        nodes << 0, 2, 2, 0, 0, 2, 2, 0,
                 0, 0, 3, 3, 0, 0, 3, 3,
                 0, 0, 0, 0, 4, 4, 4, 4;


        std::vector<int> nodeIds = myStructureBrick.NodesCreate(nodes);
        std::vector<int> ids(8);

        // right numbering
        ids = {nodeIds[0], nodeIds[1], nodeIds[2], nodeIds[3], nodeIds[4], nodeIds[5], nodeIds[6], nodeIds[7]};
        myStructureBrick.ElementCreate(itBricks, ids);

        // "wrong" numbering
        ids = {nodeIds[4], nodeIds[5], nodeIds[6], nodeIds[7], nodeIds[0], nodeIds[1], nodeIds[2], nodeIds[3]};
        myStructureBrick.ElementCreate(itBricks, ids);
    }
    catch (NuTo::MechanicsException& e)
    {
        std::cout << "Brick failed." << std::endl;
        std::cout << e.ErrorMessage() << std::endl;
        errorsOccurred = true;
    }


    if (errorsOccurred)
        throw NuTo::MechanicsException("[NodeReordering] Errors occurred!");
}


//#define TRYCATCH

int main(int argc, char* argv[])
{

    boost::filesystem::path binaryPath = std::string(argv[0]);
    binaryPath.remove_filename();

    std::string meshFile2D = binaryPath.string() + "/2D.msh";

#ifdef TRYCATCH
    try
    {
#endif
    CheckTriangle();
    CheckQuad();
    CheckTetrahedron();
    CheckBrick();
    CheckPrism();
    CheckTruss();
    CheckAPI();
    ImportFromGmsh(meshFile2D);
    NodeReordering();
#ifdef TRYCATCH
}
catch (NuTo::Exception& e)
{
    std::cout << e.ErrorMessage();
    return -1;
}
catch (...)
{
    std::cout << "Something else went wrong." << std::endl;
    return -1;
}
#endif

    std::cout << std::endl;
    std::cout << "#####################################" << std::endl;
    std::cout << "##  Congratulations! Everything    ##" << std::endl;
    std::cout << "##   went better than expected!    ##" << std::endl;
    std::cout << "#####################################" << std::endl;

    return 0;
}

