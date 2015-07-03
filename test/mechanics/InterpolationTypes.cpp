/*
 * Interpolation2D.cpp
 *
 *  Created on: 20 Mar 2015
 *      Author: ttitsche
 */

#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/interpolationtypes/InterpolationType.h"
#include "nuto/mechanics/interpolationtypes/Interpolation2DTriangle.h"


#include "nuto/mechanics/integrationtypes/IntegrationType1D2NGauss2Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType1D2NGauss3Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType1D2NGauss4Ip.h"

#include "nuto/mechanics/integrationtypes/IntegrationType2D3NGauss13Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D4NGauss4Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType3D4NGauss4Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType3D8NGauss2x2x2Ip.h"

#include "nuto/base/Exception.h"
#include <boost-1_55/boost/filesystem.hpp>

#include <map>
#include <string>
#include <iostream>

#include <cmath>

const bool PRINTRESULT = true;

//! @brief checks, whether or not the natural node coordinates match the shape functions
//! This should be true: N_i(xi_j) == 1 for i == j    and      N_j(xi_i) == 0 for i != j
void CheckShapeFunctionsAndNodePositions(NuTo::InterpolationType& rIT, int rNumNodesExpected)
{
    auto dofType = NuTo::Node::COORDINATES;
    assert(rIT.IsDof(dofType)); // coordinates exist?

    if (rIT.GetNumNodes() != rNumNodesExpected)
    {
        std::cout << rIT.GetNumNodes() << std::endl;
        std::cout << rIT.Info() << std::endl;
        throw NuTo::MechanicsException("[CheckTriangle] Wrong node number");
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
                    throw NuTo::MechanicsException("[CheckShapeFunctionsAndNodePositions] shape functions and node positions do not match (should be 1).");
            } else
            {
                // should be 0
                if (std::abs(shapeFunctions(iShapeFunctions)) > 1.e-10)
                    throw NuTo::MechanicsException("[CheckShapeFunctionsAndNodePositions] shape functions and node positions do not match (should be 0).");
            }
        }
    }
    std::cout << "[CheckShapeFunctionsAndNodePositions] OK!" << std::endl;
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

    {
        NuTo::IntegrationType1D2NGauss2Ip myIntegrationType;
        NuTo::InterpolationType myIT(NuTo::Interpolation::eShapeType::TRUSS1D);
        myIT.AddDofInterpolation(NuTo::Node::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
        CheckShapeFunctionsAndNodePositions(myIT, 3);
        myIT.UpdateIntegrationType(myIntegrationType);
        myIT.PrintNodeCoordinates();
        myIT.PrintNodeIndices();
        CheckNodeIndexing (myIT);
    }


    {
        NuTo::IntegrationType1D2NGauss3Ip myIntegrationType;
        NuTo::InterpolationType myIT(NuTo::Interpolation::eShapeType::TRUSS1D);
        myIT.AddDofInterpolation(NuTo::Node::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT3);
        CheckShapeFunctionsAndNodePositions(myIT, 4);
        myIT.UpdateIntegrationType(myIntegrationType);
        myIT.PrintNodeCoordinates();
        myIT.PrintNodeIndices();
        CheckNodeIndexing (myIT);
    }


    {
        NuTo::IntegrationType1D2NGauss4Ip myIntegrationType;
        NuTo::InterpolationType myIT(NuTo::Interpolation::eShapeType::TRUSS1D);
        myIT.AddDofInterpolation(NuTo::Node::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
        myIT.AddDofInterpolation(NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT4);
        myIT.AddDofInterpolation(NuTo::Node::NONLOCALEQSTRAIN, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
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

    NuTo::InterpolationType myIT3(NuTo::Interpolation::eShapeType::TRIANGLE2D);
    myIT3.AddDofInterpolation(NuTo::Node::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    myIT3.UpdateIntegrationType(myIntegrationType);
    CheckShapeFunctionsAndNodePositions(myIT3, 3);

    NuTo::InterpolationType myIT6(NuTo::Interpolation::eShapeType::TRIANGLE2D);
    myIT6.AddDofInterpolation(NuTo::Node::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    myIT6.UpdateIntegrationType(myIntegrationType);
    CheckShapeFunctionsAndNodePositions(myIT6, 6);

    NuTo::InterpolationType myIT10(NuTo::Interpolation::eShapeType::TRIANGLE2D);
    myIT10.AddDofInterpolation(NuTo::Node::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT3);
    myIT10.UpdateIntegrationType(myIntegrationType);
    CheckShapeFunctionsAndNodePositions(myIT10, 10);

    NuTo::InterpolationType myIT15(NuTo::Interpolation::eShapeType::TRIANGLE2D);
    myIT15.AddDofInterpolation(NuTo::Node::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT4);
    myIT15.UpdateIntegrationType(myIntegrationType);
    CheckShapeFunctionsAndNodePositions(myIT15, 15);

    NuTo::InterpolationType myIT(NuTo::Interpolation::eShapeType::TRIANGLE2D);
    myIT.AddDofInterpolation(NuTo::Node::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    myIT.AddDofInterpolation(NuTo::Node::TEMPERATURES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT4);
    myIT.AddDofInterpolation(NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT3);
    myIT.AddDofInterpolation(NuTo::Node::NONLOCALEQSTRAIN, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    myIT.UpdateIntegrationType(myIntegrationType);

    if (PRINTRESULT)
    {
        myIT.PrintNodeCoordinates();
        myIT.PrintNodeIndices();
    }

    CheckNodeIndexing(myIT);

    // check previous dofs for that specific case

    if (myIT.Get(NuTo::Node::COORDINATES).GetLocalStartIndex() != 0) // no coordinates
        throw NuTo::MechanicsException("[CheckTriangle] Wrong GetNumPreviousActiveDofs for COORDINATES");

    if (myIT.Get(NuTo::Node::DISPLACEMENTS).GetLocalStartIndex() != 0) // displacement is always first
        throw NuTo::MechanicsException("[CheckTriangle] Wrong GetNumPreviousActiveDofs for DISPLACEMENTS");

    if (myIT.Get(NuTo::Node::TEMPERATURES).GetLocalStartIndex() != 20) // 10x2 previous displacement dofs
        throw NuTo::MechanicsException("[CheckTriangle] Wrong GetNumPreviousActiveDofs for TEMPERATURE");

    if (myIT.Get(NuTo::Node::NONLOCALEQSTRAIN).GetLocalStartIndex() != 35) // + 15x1 previous temperature dofs
        throw NuTo::MechanicsException("[CheckTriangle] Wrong GetNumPreviousActiveDofs for NONLOCALEQSTRAIN");

    std::cout << "[CheckPreviousDofs] OK!" << std::endl;

}

void CheckQuad()
{

    NuTo::IntegrationType2D4NGauss4Ip myIntegrationType;

    NuTo::InterpolationType myIT4(NuTo::Interpolation::eShapeType::QUAD2D);
    myIT4.AddDofInterpolation(NuTo::Node::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    myIT4.UpdateIntegrationType(myIntegrationType);
    CheckShapeFunctionsAndNodePositions(myIT4, 4);

    NuTo::InterpolationType myIT8(NuTo::Interpolation::eShapeType::QUAD2D);
    myIT8.AddDofInterpolation(NuTo::Node::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    myIT8.UpdateIntegrationType(myIntegrationType);
    CheckShapeFunctionsAndNodePositions(myIT8, 8);

    NuTo::InterpolationType myIT9(NuTo::Interpolation::eShapeType::QUAD2D);
    myIT9.AddDofInterpolation(NuTo::Node::COORDINATES, NuTo::Interpolation::eTypeOrder::LOBATTO2);
    myIT9.UpdateIntegrationType(myIntegrationType);
    CheckShapeFunctionsAndNodePositions(myIT9, 9);

    NuTo::InterpolationType myIT16(NuTo::Interpolation::eShapeType::QUAD2D);
    myIT16.AddDofInterpolation(NuTo::Node::COORDINATES, NuTo::Interpolation::eTypeOrder::LOBATTO3);
    myIT16.UpdateIntegrationType(myIntegrationType);
    CheckShapeFunctionsAndNodePositions(myIT16, 16);

    NuTo::InterpolationType myIT25(NuTo::Interpolation::eShapeType::QUAD2D);
    myIT25.AddDofInterpolation(NuTo::Node::COORDINATES, NuTo::Interpolation::eTypeOrder::LOBATTO4);
    myIT25.UpdateIntegrationType(myIntegrationType);
    CheckShapeFunctionsAndNodePositions(myIT25, 25);

}

void CheckTetrahedron()
{

    NuTo::IntegrationType3D4NGauss4Ip myIntegrationType;

    NuTo::InterpolationType myIT4(NuTo::Interpolation::eShapeType::TETRAHEDRON3D);
    myIT4.AddDofInterpolation(NuTo::Node::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    myIT4.UpdateIntegrationType(myIntegrationType);
    CheckShapeFunctionsAndNodePositions(myIT4, 4);

    NuTo::InterpolationType myIT10(NuTo::Interpolation::eShapeType::TETRAHEDRON3D);
    myIT10.AddDofInterpolation(NuTo::Node::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    myIT10.UpdateIntegrationType(myIntegrationType);
    CheckShapeFunctionsAndNodePositions(myIT10, 10);

}

void CheckBrick()
{

    NuTo::IntegrationType3D8NGauss2x2x2Ip myIntegrationType;

    NuTo::InterpolationType myIT8(NuTo::Interpolation::eShapeType::BRICK3D);
    myIT8.AddDofInterpolation(NuTo::Node::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    myIT8.UpdateIntegrationType(myIntegrationType);
    CheckShapeFunctionsAndNodePositions(myIT8, 8);

    NuTo::InterpolationType myIT20(NuTo::Interpolation::eShapeType::BRICK3D);
    myIT20.AddDofInterpolation(NuTo::Node::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    myIT20.UpdateIntegrationType(myIntegrationType);
    CheckShapeFunctionsAndNodePositions(myIT20, 20);

}

//! @brief API of the interpolation types
void CheckAPI()
{
    NuTo::Structure myStructure(2);

    double lx = 3;
    double ly = 5;

    //create constitutive law
    int myMatLin = myStructure.ConstitutiveLawCreate("LinearElasticEngineeringStress");
    myStructure.ConstitutiveLawSetYoungsModulus(myMatLin, 10);
    myStructure.ConstitutiveLawSetPoissonsRatio(myMatLin, 0.2);

    //create section
    int mySection = myStructure.SectionCreate("Plane_Strain");
    myStructure.SectionSetThickness(mySection, 1);


    NuTo::FullVector<double, Eigen::Dynamic> nodeCoordinates(2);

    // nodes without dofs, only coordinates
    nodeCoordinates << 0, 0;
    int nodeIndex0 = myStructure.NodeCreate(nodeCoordinates);

    nodeCoordinates << lx, 0;
    int nodeIndex1 = myStructure.NodeCreate(nodeCoordinates);

    nodeCoordinates << 0, ly;
    int nodeIndex2 = myStructure.NodeCreate(nodeCoordinates);

    NuTo::FullVector<int, Eigen::Dynamic> nodeIndices(3);
    nodeIndices << nodeIndex0, nodeIndex1, nodeIndex2;

    // create interpolation type
    int myInterpolationTypeIndex = myStructure.InterpolationTypeCreate("TRIANGLE2D");

    myStructure.InterpolationTypeAdd(myInterpolationTypeIndex, NuTo::Node::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(myInterpolationTypeIndex, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    myStructure.InterpolationTypeSetIntegrationType(myInterpolationTypeIndex, NuTo::IntegrationType::IntegrationType2D3NGauss13Ip, NuTo::IpData::STATICDATA);

    int elementIndex = myStructure.ElementCreate(myInterpolationTypeIndex, nodeIndices, "ConstitutiveLawIp", "StaticData");
    myStructure.ElementSetConstitutiveLaw(elementIndex, myMatLin);

    int elementGroup = myStructure.GroupCreate("ELEMENTS");
    myStructure.GroupAddElement(elementGroup, elementIndex);

    myStructure.ElementConvertToInterpolationType(elementGroup, 0.01, 2.);

    int numIP = myStructure.ElementGetElementPtr(elementIndex)->GetNumIntegrationPoints();
    std::cout << numIP << std::endl;
    assert(numIP == 13);

    //assign constitutive law
    myStructure.ElementTotalSetSection(mySection);

    NuTo::SparseMatrixCSRVector2General<double> hessian;
    NuTo::FullVector<double, Eigen::Dynamic> internalGradient, dummy;

    myStructure.NodeBuildGlobalDofs();
    myStructure.CalculateMaximumIndependentSets();

    myStructure.BuildGlobalCoefficientMatrix0(hessian, dummy);
    myStructure.BuildGlobalGradientInternalPotentialVector(internalGradient);

    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> hessian_full(hessian);

//    hessian_full.Info(10,3, true);
//    internalGradient.Info(10,3, true);

}

//! @brief Imports a mesh file and builds the hessian and the internal gradient.
void ImportFromGmsh(std::string rMeshFile)
{
    NuTo::Structure myStructure(2);
    NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> groupIndices = myStructure.ImportFromGmsh(rMeshFile, NuTo::ElementData::CONSTITUTIVELAWIP, NuTo::IpData::NOIPDATA);

    std::cout << groupIndices.size() << std::endl;

    int interpolationType = groupIndices.GetValue(0, 1);
    myStructure.InterpolationTypeAdd(interpolationType, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::EQUIDISTANT2);

    myStructure.SetVerboseLevel(10);
    myStructure.ElementConvertToInterpolationType(groupIndices.GetValue(0, 0));

    myStructure.InterpolationTypeSetIntegrationType(interpolationType, NuTo::IntegrationType::IntegrationType2D3NGauss3Ip, NuTo::IpData::NOIPDATA);

    myStructure.InterpolationTypeInfo(0);

    myStructure.NodeBuildGlobalDofs();
    int section = myStructure.SectionCreate("PLANE_STRAIN");
    myStructure.SectionSetThickness(section, 3);

    myStructure.ElementTotalSetSection(section);

    int constitutiveLaw = myStructure.ConstitutiveLawCreate("LinearElasticEngineeringStress");
    myStructure.ConstitutiveLawSetYoungsModulus(constitutiveLaw, 42);
    myStructure.ConstitutiveLawSetPoissonsRatio(constitutiveLaw, .25);
    myStructure.ElementTotalSetConstitutiveLaw(constitutiveLaw);

    myStructure.NodeBuildGlobalDofs();
    myStructure.CalculateMaximumIndependentSets();
    NuTo::FullVector<double, Eigen::Dynamic> actDofValues, depDofValues, internalGradient, dummy;
    NuTo::SparseMatrixCSRVector2General<double> hessian0_CSR;

    myStructure.BuildGlobalGradientInternalPotentialVector(internalGradient);
    myStructure.BuildGlobalCoefficientMatrix0(hessian0_CSR, dummy);

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
    int itTruss       = myStructureTruss.InterpolationTypeCreate(NuTo::Interpolation::TRUSS1D);
    int itTriangle    = myStructureTriangle.InterpolationTypeCreate(NuTo::Interpolation::TRIANGLE2D);
    int itQuad        = myStructureQuad.InterpolationTypeCreate(NuTo::Interpolation::QUAD2D);
    int itTetrahedron = myStructureTetrahedron.InterpolationTypeCreate(NuTo::Interpolation::TETRAHEDRON3D);
    int itBricks      = myStructureBrick.InterpolationTypeCreate(NuTo::Interpolation::BRICK3D);

    myStructureTruss.InterpolationTypeAdd(itTruss, NuTo::Node::COORDINATES, NuTo::Interpolation::EQUIDISTANT1);
    myStructureTriangle.InterpolationTypeAdd(itTriangle, NuTo::Node::COORDINATES, NuTo::Interpolation::EQUIDISTANT1);
    myStructureQuad.InterpolationTypeAdd(itQuad, NuTo::Node::COORDINATES, NuTo::Interpolation::EQUIDISTANT1);
    myStructureTetrahedron.InterpolationTypeAdd(itTetrahedron, NuTo::Node::COORDINATES, NuTo::Interpolation::EQUIDISTANT1);
    myStructureBrick.InterpolationTypeAdd(itBricks, NuTo::Node::COORDINATES, NuTo::Interpolation::EQUIDISTANT1);



    NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> nodeIds;


    // **************
    // **  TRUSS
    // **************
    try
    {
        NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> nodes(1,2);
        nodes << 0, 1;

        nodeIds = myStructureTruss.NodesCreate(nodes);
        NuTo::FullVector<int, Eigen::Dynamic> ids(2);

        // right numbering
        ids << nodeIds(0), nodeIds(1);
        myStructureTruss.ElementCreate(itTruss, ids);

        // "wrong" numbering
        ids << nodeIds(1), nodeIds(0);
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

        NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> nodes(2,3);
        nodes << 0, 2, 2,
                 0, 0, 3;

        nodeIds = myStructureTriangle.NodesCreate(nodes);
        NuTo::FullVector<int, Eigen::Dynamic> ids(3);

        // right numbering
        ids << nodeIds(0), nodeIds(1), nodeIds(2);
        myStructureTriangle.ElementCreate(itTriangle, ids);

        // "wrong" numbering
        ids << nodeIds(2), nodeIds(1), nodeIds(0);
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
        NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> nodes(2,4);
        nodes << 0, 2, 2, 0,
                 0, 0, 3, 2;

        nodeIds = myStructureQuad.NodesCreate(nodes);
        NuTo::FullVector<int, Eigen::Dynamic> ids(4);

        // right numbering
        ids << nodeIds(0), nodeIds(1), nodeIds(2), nodeIds(3);
        myStructureQuad.ElementCreate(itQuad, ids);

        // "wrong" numbering
        ids << nodeIds(3), nodeIds(2), nodeIds(1), nodeIds(0);
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
        NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> nodes(3,4);
        nodes << 0, 2, 0, 0,
                 0, 0, 3, 0,
                 0, 0, 0, 4;


        nodeIds = myStructureTetrahedron.NodesCreate(nodes);
        NuTo::FullVector<int, Eigen::Dynamic> ids(4);

        // right numbering
        ids << nodeIds(0), nodeIds(1), nodeIds(2), nodeIds(3);
        myStructureTetrahedron.ElementCreate(itTetrahedron, ids);

        // "wrong" numbering
        ids << nodeIds(0), nodeIds(3), nodeIds(2), nodeIds(1);
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
        NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> nodes(3,8);
        nodes << 0, 2, 2, 0, 0, 2, 2, 0,
                 0, 0, 3, 3, 0, 0, 3, 3,
                 0, 0, 0, 0, 4, 4, 4, 4;


        nodeIds = myStructureBrick.NodesCreate(nodes);
        NuTo::FullVector<int, Eigen::Dynamic> ids(8);

        // right numbering
        ids << nodeIds(0), nodeIds(1), nodeIds(2), nodeIds(3), nodeIds(4), nodeIds(5), nodeIds(6), nodeIds(7);
        myStructureBrick.ElementCreate(itBricks, ids);

        // "wrong" numbering
        ids << nodeIds(4), nodeIds(5), nodeIds(6), nodeIds(7), nodeIds(0), nodeIds(1), nodeIds(2), nodeIds(3);
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


#define TRYCATCH

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

