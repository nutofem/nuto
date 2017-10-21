#include <math/SpatialContainer.h>
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

using namespace NuTo;

//! @brief helper struct that defines the element surface by the element ptr and the surface id
struct ElementSurface
{
    NuTo::ElementBase* mElement;
    int mSurface;
};


struct NodeCoordinate
{
    Eigen::VectorXd operator()(const NuTo::NodeBase* node) const
    {
        return node->Get(NuTo::Node::eDof::COORDINATES);
    }
};

using NodeTree = NuTo::SpatialContainer<NuTo::NodeBase*, NodeCoordinate>;


std::vector<NuTo::NodeBase*> GetIntersectingNodes(NuTo::Structure& s, int groupIdMaster, int groupIdSlave)
{
    int nMaster = s.GroupCreateNodeGroupFromElements(groupIdMaster);
    int nSlave = s.GroupCreateNodeGroupFromElements(groupIdSlave);
    int gNodesPrism = s.GroupIntersection(nMaster, nSlave);

    std::vector<NuTo::NodeBase*> v;
    v.reserve(s.GroupGetNumMembers(gNodesPrism));
    for (int nodeIndex : s.GroupGetMemberIds(gNodesPrism))
        v.push_back(s.NodeGetNodePtr(nodeIndex));

    s.GroupDelete(nMaster);
    s.GroupDelete(nSlave);
    s.GroupDelete(gNodesPrism);
    return v;
}

std::vector<NuTo::NodeBase*> GetSurfaceNodes(ElementSurface elementSurface)
{
    Eigen::VectorXi surfaceNodeIndices =
            elementSurface.mElement->GetInterpolationType().GetSurfaceNodeIndices(elementSurface.mSurface);
    auto numSurfaceNodes = surfaceNodeIndices.rows();
    std::vector<NuTo::NodeBase*> surfaceNodes(numSurfaceNodes);

    for (int iSurfaceNode = 0; iSurfaceNode < numSurfaceNodes; ++iSurfaceNode)
        surfaceNodes[iSurfaceNode] = elementSurface.mElement->GetNode(surfaceNodeIndices[iSurfaceNode]);

    return surfaceNodes;
}


bool IsSurface(ElementSurface elementSurface, const NodeTree& nodeTree)
{
    std::vector<NuTo::NodeBase*> surfaceNodes = GetSurfaceNodes(elementSurface);

    // check, if all surface nodes are in the node tree
    for (auto& surfaceNode : surfaceNodes)
    {
        Eigen::VectorXd coordinate = surfaceNode->Get(NuTo::Node::eDof::COORDINATES);
        if (not nodeTree.HasEntryAtCoordinate(coordinate, 1.e-10))
            return false;
    }
    return true;
}

int FindSurfaceId(NuTo::ElementBase* element, const NodeTree& nodeTree)
{
    const auto& it = element->GetInterpolationType();
    for (int iSurface = 0; iSurface < it.GetNumSurfaces(); ++iSurface)
        if (IsSurface({element, iSurface}, nodeTree))
            return iSurface;

    return -42; // error code.
}

std::vector<ElementSurface> GetElementSurfaceVector(NuTo::Structure& s, int elementGroup, const NodeTree& nodeTree)
{
    std::vector<ElementSurface> v;
    v.reserve(s.GroupGetNumMembers(elementGroup));
    for (int elementId : s.GroupGetMemberIds(elementGroup))
    {
        NuTo::ElementBase* e = s.ElementGetElementPtr(elementId);
        int surfaceId = FindSurfaceId(e, nodeTree);
        if (surfaceId >= 0)
            v.push_back({e, surfaceId});
        // else: element is not part of the surface
    }
    return v;
}

bool AreEqual(const std::vector<NuTo::NodeBase*>& r1, const std::vector<NuTo::NodeBase*>& r2)
{
    return std::is_permutation(r1.begin(), r1.end(), r2.begin(), r2.end());
}


std::vector<std::pair<ElementSurface, ElementSurface>> FindMatchingElements(NuTo::Structure& s, int groupIdMaster,
                                                                            int groupIdSlave)
{
    std::vector<NuTo::NodeBase*> gNodesPrism = GetIntersectingNodes(s, groupIdMaster, groupIdSlave);

    NodeTree nodeTree(gNodesPrism);
    std::vector<ElementSurface> eMaster = GetElementSurfaceVector(s, groupIdMaster, nodeTree);
    std::vector<ElementSurface> eSlave = GetElementSurfaceVector(s, groupIdSlave, nodeTree);

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
            throw NuTo::Exception(__PRETTY_FUNCTION__, "No matching Slaveate surface found.");
    }
    return pairs;
}

bool HasOnlyCoordinateInterpolation(NuTo::Structure& s, int gElement)
{
    for (auto* e : NuTo::MeshCompanion::GetElementVector(s, gElement))
    {
        auto dofs = e->GetInterpolationType().GetDofs();
        if (dofs.size() != 1 or *dofs.begin() != NuTo::Node::eDof::COORDINATES)
            return false;
    }
    return true;
}

std::vector<std::pair<ElementSurface, ElementSurface>>
FindElementPairsContainingTheNode(const NuTo::NodeBase* node,
                                  const std::vector<std::pair<ElementSurface, ElementSurface>>& surfaceElementPairs)
{
    std::vector<std::pair<ElementSurface, ElementSurface>> matchingPairs;
    for (const auto& pair : surfaceElementPairs)
    {
        NuTo::ElementBase* element = pair.first.mElement;
        for (int i = 0; i < element->GetNumNodes(); ++i)
        {
            if (element->GetNode(i) == node)
            {
                matchingPairs.push_back(pair);
                continue;
            }
        }
    }
    return matchingPairs;
}

int GetNodeCoordinatesIndex(const NuTo::NodeBase* node, const NuTo::ElementBase* element)
{
    for (int i = 0; i < element->GetNumNodes(NuTo::Node::eDof::COORDINATES); ++i)
    {
        if (element->GetNode(i, NuTo::Node::eDof::COORDINATES) == node)
            return i;
    }
    throw;
}

Eigen::VectorXd GetLocalSurfaceCoordinates(int nodeIndex, const ElementSurface& elementSurface)
{
    const auto& it = elementSurface.mElement->GetInterpolationType().Get(NuTo::Node::eDof::COORDINATES);

    // find surface parameters
    auto localNodeCoordinates = it.GetNaturalNodeCoordinates(nodeIndex);
    Eigen::VectorXd R = it.CalculateNaturalSurfaceCoordinates(Eigen::VectorXd::Zero(2), elementSurface.mSurface) -
                        localNodeCoordinates;
    Eigen::MatrixXd dRdS =
            it.CalculateDerivativeNaturalSurfaceCoordinates(Eigen::VectorXd::Zero(2), elementSurface.mSurface);

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(dRdS, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::VectorXd surfaceParameters = -svd.solve(R);
    assert((localNodeCoordinates - it.CalculateNaturalSurfaceCoordinates(surfaceParameters, elementSurface.mSurface))
                   .norm() < 1.e-10);
    return surfaceParameters;
}

Eigen::VectorXd GetLocalSurfaceCoordinates(const NuTo::NodeBase* node, const ElementSurface& elementSurface)
{
    return GetLocalSurfaceCoordinates(GetNodeCoordinatesIndex(node, elementSurface.mElement), elementSurface);
}


Eigen::VectorXd CalculateNormalAtNode(const NuTo::NodeBase* node, const ElementSurface& elementSurface)
{
    Eigen::MatrixXd nodeCoordinates = elementSurface.mElement->ExtractNodeValues(0, NuTo::Node::eDof::COORDINATES);
    const auto& it = elementSurface.mElement->GetInterpolationType().Get(NuTo::Node::eDof::COORDINATES);
    Eigen::VectorXd ipCoordsSurface = GetLocalSurfaceCoordinates(node, elementSurface);
    Eigen::VectorXd ipCoordsNatural = it.CalculateNaturalSurfaceCoordinates(ipCoordsSurface, elementSurface.mSurface);

    Eigen::MatrixXd derivativeShapeFunctionsNatural = it.DerivativeShapeFunctionsNatural(ipCoordsNatural);
    const Eigen::Matrix3d jacobian = dynamic_cast<ContinuumElement<3>*>(elementSurface.mElement)
                                             ->CalculateJacobian(derivativeShapeFunctionsNatural, nodeCoordinates);

    Eigen::MatrixXd derivativeNaturalSurfaceCoordinates = it.CalculateDerivativeNaturalSurfaceCoordinates(
            ipCoordsSurface, elementSurface.mSurface); // = [dXi / dAlpha]
    Eigen::Vector3d dXdAlpha = jacobian * derivativeNaturalSurfaceCoordinates.col(0);
    Eigen::Vector3d dXdBeta = jacobian * derivativeNaturalSurfaceCoordinates.col(1);

    Eigen::Vector3d surfaceNormalVector = dXdAlpha.cross(dXdBeta); // = || [dX / dXi] * [dXi / dAlpha] ||
    surfaceNormalVector.normalize();
    return surfaceNormalVector;
}

NuTo::NodeBase* CloneNode(NuTo::Structure& s, const NuTo::NodeBase* node)
{
    int nodeIndex = s.NodeCreate(node->Get(NuTo::Node::eDof::COORDINATES));
    return s.NodeGetNodePtr(nodeIndex);
}

NuTo::Interpolation::eTypeOrder GetCoordinateInterpolation(NuTo::Structure& s, int groupIdMaster, int)
{
    auto* firstElement = s.ElementGetElementPtr(s.GroupGetMemberIds(groupIdMaster)[0]);
    NuTo::Interpolation::eTypeOrder coordinateInterpolation =
            firstElement->GetInterpolationType().Get(NuTo::Node::eDof::COORDINATES).GetTypeOrder();
    for (int elementId : s.GroupGetMemberIds(groupIdMaster))
    {
        auto* e = s.ElementGetElementPtr(elementId);
        auto type = e->GetInterpolationType().Get(NuTo::Node::eDof::COORDINATES).GetTypeOrder();
        if (type != coordinateInterpolation)
            throw NuTo::Exception(__PRETTY_FUNCTION__,
                                  "All elements in the groups must have the same coordinate interpolation.");
    }
    return coordinateInterpolation;
}

NuTo::NodeBase* InterpolateNode(NuTo::Structure& s, NuTo::NodeBase* nodeMaster, NuTo::NodeBase* nodeSlave,
                                std::map<NuTo::NodeBase*, NuTo::NodeBase*>& masterNodeToInterpolatedNode)
{
    auto it = masterNodeToInterpolatedNode.find(nodeMaster);
    if (it != masterNodeToInterpolatedNode.end())
        return it->second;

    Eigen::Vector3d coords =
            (nodeMaster->Get(NuTo::Node::eDof::COORDINATES) + nodeSlave->Get(NuTo::Node::eDof::COORDINATES)) / 2.;
    int nodeIndex = s.NodeCreate(coords);

    NuTo::NodeBase* newNode = s.NodeGetNodePtr(nodeIndex);
    masterNodeToInterpolatedNode[nodeMaster] = newNode;
    return newNode;
}

int CreateLinearPrism(NuTo::Structure& s, ElementSurface elementSurfaceMaster, int interpolationTypeId,
                      const std::map<NuTo::NodeBase*, NuTo::NodeBase*>& clonedNodesMapping)
{
    auto surfaceNodes = GetSurfaceNodes(elementSurfaceMaster);
    // linear elements
    assert(surfaceNodes.size() == 3);
    for (int i = 0; i < 3; ++i)
    {
        NuTo::NodeBase* masterNode = surfaceNodes[i];
        NuTo::NodeBase* slaveNode = clonedNodesMapping.at(masterNode);
        surfaceNodes.push_back(slaveNode);
    }
    return s.ElementCreate(interpolationTypeId, surfaceNodes);
}

int CreateQuadraticPrism(NuTo::Structure& s, ElementSurface elementSurfaceMaster, int interpolationTypeId,
                         const std::map<NuTo::NodeBase*, NuTo::NodeBase*>& clonedNodesMapping,
                         std::map<NuTo::NodeBase*, NuTo::NodeBase*>& masterNodeToInterpolatedNode)
{
    const auto& itBase = elementSurfaceMaster.mElement->GetInterpolationType().Get(NuTo::Node::eDof::COORDINATES);

    std::vector<NuTo::NodeBase*> sortedNodes(6);

    for (int i = 0; i < itBase.GetNumSurfaceNodes(elementSurfaceMaster.mSurface); ++i)
    {
        int nodeIndex = itBase.GetSurfaceNodeIndex(elementSurfaceMaster.mSurface, i);
        Eigen::Vector2d surfaceCoords = GetLocalSurfaceCoordinates(nodeIndex, elementSurfaceMaster);

        if (surfaceCoords.isApprox(Eigen::Vector2d{0, 0}))
            sortedNodes[0] = elementSurfaceMaster.mElement->GetNode(nodeIndex);
        if (surfaceCoords.isApprox(Eigen::Vector2d{1, 0}))
            sortedNodes[1] = elementSurfaceMaster.mElement->GetNode(nodeIndex);
        if (surfaceCoords.isApprox(Eigen::Vector2d{0, 1}))
            sortedNodes[2] = elementSurfaceMaster.mElement->GetNode(nodeIndex);
        if (surfaceCoords.isApprox(Eigen::Vector2d{.5, 0}))
            sortedNodes[3] = elementSurfaceMaster.mElement->GetNode(nodeIndex);
        if (surfaceCoords.isApprox(Eigen::Vector2d{.5, .5}))
            sortedNodes[4] = elementSurfaceMaster.mElement->GetNode(nodeIndex);
        if (surfaceCoords.isApprox(Eigen::Vector2d{0, .5}))
            sortedNodes[5] = elementSurfaceMaster.mElement->GetNode(nodeIndex);
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
    nodeVector[3] = clonedNodesMapping.at(sortedNodes[0]);
    nodeVector[4] = clonedNodesMapping.at(sortedNodes[1]);
    nodeVector[5] = clonedNodesMapping.at(sortedNodes[2]);

    // lower mids
    nodeVector[6] = sortedNodes[3];
    nodeVector[9] = sortedNodes[4];
    nodeVector[7] = sortedNodes[5];

    // upper mids
    nodeVector[12] = clonedNodesMapping.at(sortedNodes[3]);
    nodeVector[14] = clonedNodesMapping.at(sortedNodes[4]);
    nodeVector[13] = clonedNodesMapping.at(sortedNodes[5]);

    // mid corners
    nodeVector[8] = InterpolateNode(s, nodeVector[0], nodeVector[3], masterNodeToInterpolatedNode);
    nodeVector[10] = InterpolateNode(s, nodeVector[1], nodeVector[4], masterNodeToInterpolatedNode);
    nodeVector[11] = InterpolateNode(s, nodeVector[2], nodeVector[5], masterNodeToInterpolatedNode);

    // mid mids
    nodeVector[15] = InterpolateNode(s, nodeVector[6], nodeVector[12], masterNodeToInterpolatedNode);
    nodeVector[16] = InterpolateNode(s, nodeVector[7], nodeVector[13], masterNodeToInterpolatedNode);
    nodeVector[17] = InterpolateNode(s, nodeVector[9], nodeVector[14], masterNodeToInterpolatedNode);

    return s.ElementCreate(interpolationTypeId, nodeVector);
}


std::pair<int, int> NuTo::MeshCompanion::ElementPrismsCreate(NuTo::Structure& s, int groupIdMaster, int groupIdSlave,
                                                             double thickness)
{
    Timer timer(__FUNCTION__, s.GetShowTime(), s.GetLogger());

    if (not HasOnlyCoordinateInterpolation(s, groupIdMaster) or not HasOnlyCoordinateInterpolation(s, groupIdSlave))
        throw NuTo::Exception(__PRETTY_FUNCTION__, "Elements must only have COORDINATES interpolation.");

    auto pairs = FindMatchingElements(s, groupIdMaster, groupIdSlave);

    // UPDATE NODES:
    std::vector<NuTo::NodeBase*> gNodesPrism = GetIntersectingNodes(s, groupIdMaster, groupIdSlave);
    std::map<NuTo::NodeBase*, NuTo::NodeBase*> clonedNodesMapping;
    std::map<NuTo::NodeBase*, Eigen::Vector3d> normals;

    for (NuTo::NodeBase* node : gNodesPrism)
    {
        auto elementPairsContainingTheNode = FindElementPairsContainingTheNode(node, pairs);
        normals[node] = CalculateNormalAtNode(node, elementPairsContainingTheNode[0].first);
    }

    for (NuTo::NodeBase* node : gNodesPrism)
    {
        NuTo::NodeBase* clonedNode = CloneNode(s, node);
        Eigen::Vector3d nodeCoordinates = node->Get(NuTo::Node::eDof::COORDINATES);
        const Eigen::Vector3d& normal = normals[node];

        node->Set(NuTo::Node::eDof::COORDINATES, nodeCoordinates - thickness * 0.5 * normal);
        clonedNode->Set(NuTo::Node::eDof::COORDINATES, nodeCoordinates + thickness * 0.5 * normal);

        clonedNodesMapping[node] = clonedNode;

        for (int slaveElementId : s.GroupGetMemberIds(groupIdSlave))
            s.ElementGetElementPtr(slaveElementId)->ExchangeNodePtr(node, clonedNode);
    }

    std::map<NuTo::NodeBase*, NuTo::NodeBase*> masterNodeToInterpolatedNode;
    int gPrism = s.GroupCreate(NuTo::eGroupId::Elements);


    NuTo::Interpolation::eTypeOrder coordinateInterpolation =
            GetCoordinateInterpolation(s, groupIdMaster, groupIdSlave);
    int it = s.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::PRISM3D);
    s.InterpolationTypeAdd(it, NuTo::Node::eDof::COORDINATES, coordinateInterpolation);

    // CREATE PRISMS
    for (auto& pair : pairs)
    {
        if (coordinateInterpolation == NuTo::Interpolation::eTypeOrder::EQUIDISTANT1)
        {
            s.GroupAddElement(gPrism, CreateLinearPrism(s, pair.first, it, clonedNodesMapping));
        }
        else if (coordinateInterpolation == NuTo::Interpolation::eTypeOrder::EQUIDISTANT2)
        {
            s.GroupAddElement(
                    gPrism, CreateQuadraticPrism(s, pair.first, it, clonedNodesMapping, masterNodeToInterpolatedNode));
        }
        else
        {
            throw NuTo::Exception(__PRETTY_FUNCTION__,
                                  "Only implemented for EQUIDISTANT1 and EQUIDISTANT2 coordinate interpolation");
        }
    }
    return std::make_pair(gPrism, it);
}
