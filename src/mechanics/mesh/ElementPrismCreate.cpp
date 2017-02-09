#include <math/SpacialContainer.h>
#include "mechanics/mesh/MeshCompanion.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/elements/ElementBase.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/interpolationtypes/InterpolationType.h"
#include "mechanics/interpolationtypes/InterpolationBase.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "mechanics/groups/GroupEnum.h"
#include "mechanics/elements/ContinuumElement.h"
#include "base/Timer.h"

//! @brief helper struct that defines the element surface by the element ptr and the surface id
struct ElementSurface
{
    NuTo::ElementBase* mElement;
    int mSurface;
};


struct NodeCoordinate
{
    Eigen::VectorXd operator () (const NuTo::NodeBase* rNode)
    {
        return rNode->Get(NuTo::Node::eDof::COORDINATES);
    }
};

using NodeTree = NuTo::SpacialContainer<NuTo::NodeBase*, NodeCoordinate>;


std::vector<NuTo::NodeBase*> GetIntersectingNodes(NuTo::Structure &rS, int rGroupMaster, int rGroupSlave)
{
    int nMaster = rS.GroupCreateNodeGroupFromElements(rGroupMaster);
    int nSlave = rS.GroupCreateNodeGroupFromElements(rGroupSlave);
    int gNodesPrism = rS.GroupIntersection(nMaster, nSlave);

    std::vector<NuTo::NodeBase*> v;
    v.reserve(rS.GroupGetNumMembers(gNodesPrism));
    for (int nodeId : rS.GroupGetMemberIds(gNodesPrism))
        v.push_back(rS.NodeGetNodePtr(nodeId));

    rS.GroupDelete(nMaster);
    rS.GroupDelete(nSlave);
    rS.GroupDelete(gNodesPrism);
    return v;
}

std::vector<NuTo::NodeBase*> GetSurfaceNodes(ElementSurface rElementSurface)
{
    Eigen::VectorXi surfaceNodeIndices = rElementSurface.mElement->GetInterpolationType().GetSurfaceNodeIndices(rElementSurface.mSurface);
    auto numSurfaceNodes = surfaceNodeIndices.rows();
    std::vector<NuTo::NodeBase *> surfaceNodes(numSurfaceNodes);

    for (int iSurfaceNode = 0; iSurfaceNode < numSurfaceNodes; ++iSurfaceNode)
        surfaceNodes[iSurfaceNode] = rElementSurface.mElement->GetNode(surfaceNodeIndices(iSurfaceNode, 0));

    return surfaceNodes;
}



bool IsSurface(ElementSurface rElementSurface, const NodeTree& rNodeTree)
{
    Eigen::VectorXi surfaceNodeIndices = rElementSurface.mElement->GetInterpolationType().GetSurfaceNodeIndices(rElementSurface.mSurface);

    std::vector<NuTo::NodeBase *> surfaceNodes = GetSurfaceNodes(rElementSurface);

    //check, if all surface nodes are in the node tree
    for (unsigned int iNode = 0; iNode < surfaceNodes.size(); iNode++)
    {
        Eigen::VectorXd coordinate = surfaceNodes[iNode]->Get(NuTo::Node::eDof::COORDINATES);
        if (not rNodeTree.HasEntryAtCoordinate(coordinate, 1.e-10))
            return false;
    }
    return true;
}

int FindSurfaceId(NuTo::ElementBase* rElement, const NodeTree& rNodeTree)
{
    const auto& it = rElement->GetInterpolationType();
    for (int iSurface = 0; iSurface < it.GetNumSurfaces(); ++iSurface)
        if (IsSurface({rElement, iSurface}, rNodeTree))
            return iSurface;

    return -42; // error code.
}

std::vector<ElementSurface> GetElementSurfaceVector(NuTo::Structure& rS, int rG, const NodeTree& rNodeTree)
{
    std::vector<ElementSurface> v;
    v.reserve(rS.GroupGetNumMembers(rG));
    for (int elementId : rS.GroupGetMemberIds(rG))
    {
        NuTo::ElementBase* e = rS.ElementGetElementPtr(elementId);
        int surfaceId = FindSurfaceId(e, rNodeTree);
        if (surfaceId >= 0)
            v.push_back({e, surfaceId});
        // else: element is not part of the surface
    }
    return v;
}

bool AreEqual(const std::vector<NuTo::NodeBase*>& r1, const std::vector<NuTo::NodeBase*>& r2)
{
    for (const auto& n1 : r1)
    {
        bool n1IsInr2 = false;
        for (const auto& n2 : r2)
        {
            if (n1 == n2)
                n1IsInr2 = true;
        }
        if (not n1IsInr2)
            return false;
    }
    return true;
}


std::vector<std::pair<ElementSurface, ElementSurface>> FindMatchingElements(
    NuTo::Structure& rS, int rGroupMaster, int rGroupSlave)
{
    std::vector<NuTo::NodeBase*> gNodesPrism = GetIntersectingNodes(rS, rGroupMaster, rGroupSlave);

    NodeTree nodeTree(gNodesPrism);
    std::vector<ElementSurface> eMaster = GetElementSurfaceVector(rS, rGroupMaster, nodeTree);
    std::vector<ElementSurface> eSlave = GetElementSurfaceVector(rS, rGroupSlave, nodeTree);

    std::vector<std::pair<ElementSurface, ElementSurface>> pairs;
    pairs.reserve(eMaster.size());

    for (auto& esMaster : eMaster)
    {
        bool surfaceFound = false;
        std::vector<NuTo::NodeBase*> surfaceNodesMaster = GetSurfaceNodes(esMaster);
        for (auto& esSlave : eSlave)
        {
            std::vector<NuTo::NodeBase*> surfaceNodesSlave = GetSurfaceNodes(esSlave);
            if (AreEqual(surfaceNodesSlave, surfaceNodesMaster))
            {
                surfaceFound = true;
                pairs.push_back(std::make_pair(esMaster, esSlave));
            }
        }
        if (not surfaceFound)
            throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "No matching Slaveate surface found.");
    }
    return pairs;
}

bool HasOnlyCoordinateInterpolation(NuTo::Structure& rS, int gElement)
{
    for (auto* e : NuTo::MeshCompanion::GetElementVector(rS, gElement))
    {
        auto dofs = e->GetInterpolationType().GetDofs();
        if (dofs.size() != 1 or *dofs.begin() != NuTo::Node::eDof::COORDINATES)
            return false;
    }
    return true;
}

std::vector<std::pair<ElementSurface, ElementSurface>> FindElementPairsContainingTheNode(
    const NuTo::NodeBase* rNode, const std::vector<std::pair<ElementSurface, ElementSurface>>& rPairs)
{
    std::vector<std::pair<ElementSurface, ElementSurface>> matchingPairs;
    for (const auto& pair : rPairs)
    {
        NuTo::ElementBase* element = pair.first.mElement;
        for (int i = 0; i < element->GetNumNodes(); ++i)
        {
            if (element->GetNode(i) == rNode)
            {
                matchingPairs.push_back(pair);
                continue;
            }
        }
    }
    return matchingPairs;
}

int GetNodeCoordinatesIndex(const NuTo::NodeBase* rNode, const NuTo::ElementBase* rElement)
{
    for (int i = 0; i < rElement->GetNumNodes(NuTo::Node::eDof::COORDINATES); ++i)
    {
        if (rElement->GetNode(i, NuTo::Node::eDof::COORDINATES) == rNode)
            return i;
    }
    throw;
}

Eigen::VectorXd GetLocalSurfaceCoordinates(int rNodeIndex, const ElementSurface& rElementSurface)
{
    const auto& it = rElementSurface.mElement->GetInterpolationType().Get(NuTo::Node::eDof::COORDINATES);

    // find surface parameters
    auto localNodeCoordinates = it.GetNaturalNodeCoordinates(rNodeIndex);
    Eigen::VectorXd R = it.CalculateNaturalSurfaceCoordinates(Eigen::VectorXd::Zero(2), rElementSurface.mSurface) - localNodeCoordinates;
    Eigen::MatrixXd dRdS = it.CalculateDerivativeNaturalSurfaceCoordinates(Eigen::VectorXd::Zero(2), rElementSurface.mSurface);

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(dRdS, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::VectorXd surfaceParameters = -svd.solve(R);
    assert((localNodeCoordinates - it.CalculateNaturalSurfaceCoordinates(surfaceParameters, rElementSurface.mSurface)).norm() < 1.e-10);
    return surfaceParameters;
}

Eigen::VectorXd GetLocalSurfaceCoordinates(const NuTo::NodeBase* rNode, const ElementSurface& rElementSurface)
{
    return GetLocalSurfaceCoordinates(GetNodeCoordinatesIndex(rNode, rElementSurface.mElement), rElementSurface);
}




Eigen::VectorXd CalculateNormalAtNode(const NuTo::NodeBase* rNode, const ElementSurface& rElementSurface)
{
    Eigen::MatrixXd nodeCoordinates = rElementSurface.mElement->ExtractNodeValues(0, NuTo::Node::eDof::COORDINATES);
    const auto& it = rElementSurface.mElement->GetInterpolationType().Get(NuTo::Node::eDof::COORDINATES);
    Eigen::VectorXd ipCoordsSurface = GetLocalSurfaceCoordinates(rNode, rElementSurface);
    Eigen::VectorXd ipCoordsNatural = it.CalculateNaturalSurfaceCoordinates(ipCoordsSurface, rElementSurface.mSurface);

    Eigen::MatrixXd derivativeShapeFunctionsNatural = it.CalculateDerivativeShapeFunctionsNatural(ipCoordsNatural);
    const Eigen::Matrix3d jacobian = rElementSurface.mElement->AsContinuumElement3D().CalculateJacobian(derivativeShapeFunctionsNatural, nodeCoordinates);

    Eigen::MatrixXd derivativeNaturalSurfaceCoordinates = it.CalculateDerivativeNaturalSurfaceCoordinates(ipCoordsSurface, rElementSurface.mSurface); // = [dXi / dAlpha]
    Eigen::Vector3d dXdAlpha = jacobian * derivativeNaturalSurfaceCoordinates.col(0);
    Eigen::Vector3d dXdBeta  = jacobian * derivativeNaturalSurfaceCoordinates.col(1);

    Eigen::Vector3d surfaceNormalVector = dXdAlpha.cross(dXdBeta); // = || [dX / dXi] * [dXi / dAlpha] ||
    surfaceNormalVector.normalize();
    return surfaceNormalVector;
}

NuTo::NodeBase* CloneNode(NuTo::Structure& rS, const NuTo::NodeBase* rNode)
{
    int nodeId = rS.NodeCreate(rNode->Get(NuTo::Node::eDof::COORDINATES));
    return rS.NodeGetNodePtr(nodeId);
}

NuTo::Interpolation::eTypeOrder GetCoordinateInterpolation(NuTo::Structure& rS, int rGroupMaster, int rGroupSlave)
{
    auto* firstElement = rS.ElementGetElementPtr(rS.GroupGetMemberIds(rGroupMaster)[0]);
    NuTo::Interpolation::eTypeOrder coordinateInterpolation = firstElement->GetInterpolationType().Get(NuTo::Node::eDof::COORDINATES).GetTypeOrder();
    for (int elementId : rS.GroupGetMemberIds(rGroupMaster))
    {
        auto* e = rS.ElementGetElementPtr(elementId);
        auto type = e->GetInterpolationType().Get(NuTo::Node::eDof::COORDINATES).GetTypeOrder();
        if (type != coordinateInterpolation)
            throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "All elements in the groups must have the same coordinate interpolation.");
    }
    return coordinateInterpolation;
}

NuTo::NodeBase* InterpolateNode(NuTo::Structure& rS, NuTo::NodeBase* rMaster, NuTo::NodeBase* rSlave,
                                std::map<NuTo::NodeBase*, NuTo::NodeBase*>& rMasterNodeToInterpolatedNode)
{
    auto it = rMasterNodeToInterpolatedNode.find(rMaster);
    if (it != rMasterNodeToInterpolatedNode.end())
        return it->second;

    Eigen::Vector3d coords = (rMaster->Get(NuTo::Node::eDof::COORDINATES) + rSlave->Get(NuTo::Node::eDof::COORDINATES)) / 2.;
    int nodeId = rS.NodeCreate(coords);

    NuTo::NodeBase* newNode = rS.NodeGetNodePtr(nodeId);
    rMasterNodeToInterpolatedNode[rMaster] = newNode;
    return newNode;
}

int CreateLinearPrism(NuTo::Structure& rS,
                      ElementSurface rElementSurfaceMaster,
                      int rInterpolationTypeId,
                      const std::map<NuTo::NodeBase*, NuTo::NodeBase*>& rClonedNodesMapping)
{
    auto surfaceNodes = GetSurfaceNodes(rElementSurfaceMaster);
    // linear elements
    assert(surfaceNodes.size() == 3);
    for (int i = 0; i < 3; ++i)
    {
        NuTo::NodeBase* masterNode = surfaceNodes[i];
        NuTo::NodeBase* slaveNode = rClonedNodesMapping.at(masterNode);
        surfaceNodes.push_back(slaveNode);
    }
    return rS.ElementCreate(rInterpolationTypeId, surfaceNodes);
}

int CreateQuadraticPrism(NuTo::Structure& rS,
                         ElementSurface rElementSurfaceMaster,
                         int rInterpolationTypeId,
                         const std::map<NuTo::NodeBase*, NuTo::NodeBase*>& rClonedNodesMapping,
                         std::map<NuTo::NodeBase*, NuTo::NodeBase*>& rMasterNodeToInterpolatedNode)
{
    const auto& itBase = rElementSurfaceMaster.mElement->GetInterpolationType().Get(NuTo::Node::eDof::COORDINATES);

    std::vector<NuTo::NodeBase*> sortedNodes(6);

    for (int i = 0; i < itBase.GetNumSurfaceNodes(rElementSurfaceMaster.mSurface); ++i)
    {
        int nodeIndex = itBase.GetSurfaceNodeIndex(rElementSurfaceMaster.mSurface, i);
        Eigen::Vector2d surfaceCoords = GetLocalSurfaceCoordinates(nodeIndex, rElementSurfaceMaster);

        if (surfaceCoords.isApprox(Eigen::Vector2d{0, 0}))
            sortedNodes[0] = rElementSurfaceMaster.mElement->GetNode(nodeIndex);
        if (surfaceCoords.isApprox(Eigen::Vector2d{1, 0}))
            sortedNodes[1] = rElementSurfaceMaster.mElement->GetNode(nodeIndex);
        if (surfaceCoords.isApprox(Eigen::Vector2d{0, 1}))
            sortedNodes[2] = rElementSurfaceMaster.mElement->GetNode(nodeIndex);
        if (surfaceCoords.isApprox(Eigen::Vector2d{.5, 0}))
            sortedNodes[3] = rElementSurfaceMaster.mElement->GetNode(nodeIndex);
        if (surfaceCoords.isApprox(Eigen::Vector2d{.5, .5}))
            sortedNodes[4] = rElementSurfaceMaster.mElement->GetNode(nodeIndex);
        if (surfaceCoords.isApprox(Eigen::Vector2d{0, .5}))
            sortedNodes[5] = rElementSurfaceMaster.mElement->GetNode(nodeIndex);
    }

    /*
                             Prism18:

                                   w
                                   ^
                                   |
                                   3
                                 ,/|`\
                               12  |  13
                             ,/    |    `\
                            4------14-----5
                            |      8      |
                            |    ,/|`\    |
                            |  15  |  16  |
                            |,/    |    `\|
                           ,10-----17-----11
                         ,/ |      0      | `\
                        u   |    ,/ `\    |   v
                            |  ,6     `7  |
                            |,/         `\|
                            1------9------2

     */

    std::vector<NuTo::NodeBase*> nodeVector(18);
    // lower corners
    nodeVector[0] = sortedNodes[0];
    nodeVector[1] = sortedNodes[1];
    nodeVector[2] = sortedNodes[2];

    // upper corners
    nodeVector[3] = rClonedNodesMapping.at(sortedNodes[0]);
    nodeVector[4] = rClonedNodesMapping.at(sortedNodes[1]);
    nodeVector[5] = rClonedNodesMapping.at(sortedNodes[2]);

    // lower mids
    nodeVector[6] = sortedNodes[3];
    nodeVector[9] = sortedNodes[4];
    nodeVector[7] = sortedNodes[5];

    // upper mids
    nodeVector[12] = rClonedNodesMapping.at(sortedNodes[3]);
    nodeVector[14] = rClonedNodesMapping.at(sortedNodes[4]);
    nodeVector[13] = rClonedNodesMapping.at(sortedNodes[5]);

    // mid corners
    nodeVector[8] = InterpolateNode(rS, nodeVector[0], nodeVector[3], rMasterNodeToInterpolatedNode);
    nodeVector[10] = InterpolateNode(rS, nodeVector[1], nodeVector[4], rMasterNodeToInterpolatedNode);
    nodeVector[11] = InterpolateNode(rS, nodeVector[2], nodeVector[5], rMasterNodeToInterpolatedNode);

    // mid mids
    nodeVector[15] = InterpolateNode(rS, nodeVector[6], nodeVector[12], rMasterNodeToInterpolatedNode);
    nodeVector[16] = InterpolateNode(rS, nodeVector[7], nodeVector[13], rMasterNodeToInterpolatedNode);
    nodeVector[17] = InterpolateNode(rS, nodeVector[9], nodeVector[14], rMasterNodeToInterpolatedNode);

    return rS.ElementCreate(rInterpolationTypeId, nodeVector);
}


std::pair<int, int> NuTo::MeshCompanion::ElementPrismsCreate(NuTo::Structure& rS, int rGroupMaster, int rGroupSlave, double rThickness)
{
    Timer timer(__FUNCTION__, rS.GetShowTime(), rS.GetLogger());

    if (not HasOnlyCoordinateInterpolation(rS, rGroupMaster) or not HasOnlyCoordinateInterpolation(rS, rGroupSlave))
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Elements must only have COORDINATES interpolation.");

    auto pairs = FindMatchingElements(rS, rGroupMaster, rGroupSlave);

    // UPDATE NODES:
    std::vector<NuTo::NodeBase*> gNodesPrism = GetIntersectingNodes(rS, rGroupMaster, rGroupSlave);
    std::map<NuTo::NodeBase*, NuTo::NodeBase*> clonedNodesMapping;
    std::map<NuTo::NodeBase*, Eigen::Vector3d> normals;

    for (NuTo::NodeBase* node : gNodesPrism)
    {
        auto elementPairsContainingTheNode = FindElementPairsContainingTheNode(node, pairs);
        normals[node] = CalculateNormalAtNode(node, elementPairsContainingTheNode[0].first);
    }

    for (NuTo::NodeBase* node : gNodesPrism)
    {
        NuTo::NodeBase* clonedNode = CloneNode(rS, node);
        Eigen::Vector3d nodeCoordinates = node->Get(NuTo::Node::eDof::COORDINATES);
        const Eigen::Vector3d& normal = normals[node];

        node->Set(NuTo::Node::eDof::COORDINATES, nodeCoordinates - rThickness * 0.5 * normal);
        clonedNode->Set(NuTo::Node::eDof::COORDINATES, nodeCoordinates + rThickness * 0.5 * normal);

        clonedNodesMapping[node] = clonedNode;

        for (int slaveElementId : rS.GroupGetMemberIds(rGroupSlave))
            rS.ElementGetElementPtr(slaveElementId)->ExchangeNodePtr(node, clonedNode);
    }

    std::map<NuTo::NodeBase*, NuTo::NodeBase*> masterNodeToInterpolatedNode;
    int gPrism = rS.GroupCreate(NuTo::eGroupId::Elements);


    NuTo::Interpolation::eTypeOrder coordinateInterpolation = GetCoordinateInterpolation(rS, rGroupMaster, rGroupSlave);
    int it = rS.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::PRISM3D);
    rS.InterpolationTypeAdd(it, NuTo::Node::eDof::COORDINATES, coordinateInterpolation);

    // CREATE PRISMS
    for (auto& pair : pairs)
    {
        if (coordinateInterpolation == NuTo::Interpolation::eTypeOrder::EQUIDISTANT1)
        {
            rS.GroupAddElement(gPrism, CreateLinearPrism(rS, pair.first, it, clonedNodesMapping));
        }
        else if (coordinateInterpolation == NuTo::Interpolation::eTypeOrder::EQUIDISTANT2)
        {
            rS.GroupAddElement(gPrism, CreateQuadraticPrism(rS, pair.first, it, clonedNodesMapping, masterNodeToInterpolatedNode));
        }
        else
        {
            throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Only implemented for EQUIDISTANT1 and EQUIDISTANT2 coordinate interpolation");
        }
    }
    return std::make_pair(gPrism, it);
}
