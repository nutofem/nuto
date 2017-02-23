#include "BoostUnitTest.h"

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

//! @brief checks if the shape functions sum up to 1 (at all integration points)
void CheckPartitionOfUnity(const NuTo::InterpolationType &rIT, const NuTo::Node::eDof &dofType)
{
    for (int iIP = 0; iIP < rIT.GetCurrentIntegrationType().GetNumIntegrationPoints(); ++iIP)
    {
        Eigen::VectorXd ip = rIT.GetCurrentIntegrationType().GetLocalIntegrationPointCoordinates(iIP);
        BOOST_CHECK_CLOSE(rIT.Get(dofType).CalculateShapeFunctions(ip).sum(),  1., 1.e-6);
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
        BOOST_CHECK_SMALL((B - B_CDF).cwiseAbs().maxCoeff(), 1.e-4);
    }
}

//! @brief checks, whether or not the natural node coordinates match the shape functions
//! This should be true: N_i(xi_j) == 1 for i == j    and      N_j(xi_i) == 0 for i != j
void CheckShapeFunctionsAndNodePositions(NuTo::InterpolationType& rIT, int rNumNodesExpected)
{
    auto dofType = NuTo::Node::eDof::COORDINATES;
    BOOST_CHECK(rIT.IsDof(dofType));
    BOOST_CHECK_EQUAL(rIT.GetNumNodes(), rNumNodesExpected);

    int numNodes = rIT.GetNumNodes();
    for (int iNode = 0; iNode < numNodes; ++iNode)
    {
        auto naturalNodeCoordinate = rIT.Get(dofType).CalculateNaturalNodeCoordinates(iNode);
        auto shapeFunctions = rIT.Get(dofType).CalculateShapeFunctions(naturalNodeCoordinate);
        BOOST_CHECK_EQUAL(shapeFunctions.rows(), numNodes);

        for (int iShapeFunctions = 0; iShapeFunctions < shapeFunctions.rows(); ++iShapeFunctions)
        {
            if (iShapeFunctions == iNode)
                BOOST_CHECK_CLOSE(shapeFunctions(iShapeFunctions), 1, 1.e-8);
            else
                BOOST_CHECK_SMALL(shapeFunctions(iShapeFunctions), 1.e-8);
        }
    }
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

            for (int iDim = 0; iDim < coordinatesGlobal.rows(); ++iDim)
                BOOST_CHECK_CLOSE(coordinatesGlobal(iDim), coordinatesLocal(iDim), 1.e-8);
        }
    }
}

BOOST_AUTO_TEST_CASE(InterpolationTruss)
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
}

BOOST_AUTO_TEST_CASE(InterpolationTriangle)
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

    CheckNodeIndexing(myIT);
}

BOOST_AUTO_TEST_CASE(InterpolationQuad)
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

BOOST_AUTO_TEST_CASE(InterpolationTetrahedron)
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

BOOST_AUTO_TEST_CASE(InterpolationBrick)
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

BOOST_AUTO_TEST_CASE(InterpolationPrism)
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


BOOST_AUTO_TEST_CASE(InterpolationNodeReorderingTruss)
{

    NuTo::Structure myStructureTruss(1);
    int itTruss = myStructureTruss.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TRUSS1D);
    myStructureTruss.InterpolationTypeAdd(itTruss,
                                          NuTo::Node::eDof::COORDINATES,
                                          NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);

    Eigen::MatrixXd nodes(1, 2);
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
BOOST_AUTO_TEST_CASE(InterpolationNodeReorderingTriangle)
{
    NuTo::Structure myStructureTriangle(2);
    int itTriangle = myStructureTriangle.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TRIANGLE2D);
    myStructureTriangle.InterpolationTypeAdd(itTriangle,
                                             NuTo::Node::eDof::COORDINATES,
                                             NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);

    Eigen::MatrixXd nodes(2, 3);
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

BOOST_AUTO_TEST_CASE(InterpolationNodeReorderingQuad)
{
    NuTo::Structure myStructureQuad(2);
    int itQuad = myStructureQuad.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::QUAD2D);
    myStructureQuad.InterpolationTypeAdd(itQuad,
                                         NuTo::Node::eDof::COORDINATES,
                                         NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);

    Eigen::MatrixXd nodes(2, 4);
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
BOOST_AUTO_TEST_CASE(InterpolationNodeReorderingTetrahedron)
{
    NuTo::Structure myStructureTetrahedron(3);
    int itTetrahedron = myStructureTetrahedron.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TETRAHEDRON3D);
    myStructureTetrahedron.InterpolationTypeAdd(itTetrahedron,
                                                NuTo::Node::eDof::COORDINATES,
                                                NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);

    Eigen::MatrixXd nodes(3, 4);
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
BOOST_AUTO_TEST_CASE(InterpolationNodeReorderingBrick)
{
    NuTo::Structure myStructureBrick(3);
    int itBricks = myStructureBrick.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::BRICK3D);
    myStructureBrick.InterpolationTypeAdd(itBricks,
                                          NuTo::Node::eDof::COORDINATES,
                                          NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);

    Eigen::MatrixXd nodes(3, 8);
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