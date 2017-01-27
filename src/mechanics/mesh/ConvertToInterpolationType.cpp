#include <mechanics/nodes/NodeEnum.h>
#include "mechanics/mesh/MeshCompanion.h"
#include "mechanics/structures/unstructured/Structure.h"

#include "mechanics/elements/ElementBase.h"
#include "mechanics/nodes/NodeBase.h"
#include "base/Timer.h"
#include "mechanics/interpolationtypes/InterpolationType.h"
#include "mechanics/interpolationtypes/InterpolationBase.h"
#include <map>


double GetSmallestElementSize(NuTo::Structure& rS, int rGroupNumberElements)
{
    // calculate and store the 'size' of each element
    // = volume in 3D
    // = area in 2D
    // = length in 1D
    double smallestElementSize = 1e42;
    for (int elementId : rS.GroupGetMemberIds(rGroupNumberElements))
    {
        NuTo::ElementBase *element = rS.ElementGetElementPtr(elementId);
        Eigen::VectorXd sizeForEachIntegrationPoint = element->GetIntegrationPointVolume();
        smallestElementSize = std::min(smallestElementSize, sizeForEachIntegrationPoint.sum());
    }
    return smallestElementSize;
}

void NuTo::MeshCompanion::ElementTotalConvertToInterpolationType(Structure& rS)
{
    int groupNumber = rS.GroupGetElementsTotal();
    ElementConvertToInterpolationType(rS, groupNumber);
    rS.GroupDelete(groupNumber);
}

void NuTo::MeshCompanion::ElementConvertToInterpolationType(Structure& rS, int rGroupNumberElements)
{
    double sizeMin = GetSmallestElementSize(rS, rGroupNumberElements);
    double lengthMin = std::pow(sizeMin, 1. / rS.GetDimension());
    double mergeDist = lengthMin / 1000.;
    ElementConvertToInterpolationType(rS, rGroupNumberElements, mergeDist);
}


void NuTo::MeshCompanion::ElementTotalConvertToInterpolationType(Structure& rS, double rNodeDistanceMerge)
{
    int groupNumber = rS.GroupGetElementsTotal();
    ElementConvertToInterpolationType(rS, groupNumber, rNodeDistanceMerge);
    rS.GroupDelete(groupNumber);
}

struct NodeElementPair
{
    NuTo::NodeBase* node;
    NuTo::ElementBase* element;
};

struct OriginalNodeElement
{
    NuTo::NodeBase* node = nullptr;
    std::vector<NuTo::ElementBase*> elements;
};

struct SameCoordinateNodes
{
    OriginalNodeElement original;
    std::vector<NodeElementPair> newEntries;
};

struct CoordinateCompare
{
    CoordinateCompare(double rTolerance) : mTolerance(rTolerance) {}
    bool operator() (const Eigen::VectorXd& r1, const Eigen::VectorXd& r2)
    {
        for (int i = 0; i < r1.rows(); ++i)
            if (std::abs(r1[i] - r2[i]) > mTolerance)
                return (r1[i] > r2[i]);
        return false;
    }
private:
    const double mTolerance;
};

using NodeMap = std::map<Eigen::VectorXd, SameCoordinateNodes, CoordinateCompare>;

void AddExistingNodesToNodesMap(NuTo::Structure& rS, NodeMap& rNodeMap, std::vector<NuTo::ElementBase*> rElements)
{
    // add existing nodes to the node map as "orignalNode"
    for (NuTo::ElementBase* element : rElements)
    {
        // loop through nodes with coordinates
        for (int iNode = 0; iNode < element->GetNumNodes(NuTo::Node::eDof::COORDINATES); ++iNode)
        {
            NuTo::NodeBase* node = element->GetNode(iNode, NuTo::Node::eDof::COORDINATES);
            Eigen::VectorXd nodeCoordinates = node->Get(NuTo::Node::eDof::COORDINATES);
            assert(nodeCoordinates.rows() == rS.GetDimension());
            rNodeMap[nodeCoordinates].original.node = node;
            rNodeMap[nodeCoordinates].original.elements.push_back(element);
        }
    }
}

void CreateAndAddNewNodesToNodeMap(NuTo::Structure& rS, NodeMap& rNodeMap, std::vector<NuTo::ElementBase*> rElements)
{
    for (NuTo::ElementBase* element : rElements)
    {
        const NuTo::InterpolationType& interpolationType = element->GetInterpolationType();
        element->ResizeNodes(interpolationType.GetNumNodes());
        for (int iNode = 0; iNode < element->GetNumNodes(); ++iNode)
        {
            Eigen::VectorXd naturalNodeCoordinates = interpolationType.GetNaturalNodeCoordinates(iNode);
            Eigen::VectorXd globalNodeCoordinates = element->InterpolateDofGlobal(naturalNodeCoordinates, NuTo::Node::eDof::COORDINATES);

            const std::set<NuTo::Node::eDof>& nodeDofs = interpolationType.GetNodeDofs(iNode);
            NuTo::NodeBase* node = rS.NodePtrCreate(nodeDofs, globalNodeCoordinates);

            rNodeMap[globalNodeCoordinates].newEntries.push_back({node, element});
            element->SetNode(iNode, node);
        }
    }
}

void NuTo::MeshCompanion::ElementConvertToInterpolationType(Structure& rS, int rGroupNumberElements, double rNodeDistanceMerge)
{
    Timer timer(__FUNCTION__, rS.GetShowTime(), rS.GetLogger());

    std::vector<ElementBase*> elements = GetElementVector(rS, rGroupNumberElements);


    CoordinateCompare compare(rNodeDistanceMerge);
    NodeMap nodeMap(compare);

    AddExistingNodesToNodesMap(rS, nodeMap, elements);

    CreateAndAddNewNodesToNodeMap(rS, nodeMap, elements);

    std::map<const NodeBase*, int> nodeToId;
    int gAllNodes = rS.GroupGetNodesTotal();
    for (int nodeId : rS.GroupGetMemberIds(gAllNodes))
        nodeToId[rS.NodeGetNodePtr(nodeId)] = nodeId;

    // clean up duplicate nodeMap
    for (auto& coordinateNodePair : nodeMap)
    {
        SameCoordinateNodes& sameCoordinateNodes = coordinateNodePair.second;
        const auto& dofTypes = sameCoordinateNodes.newEntries[0].node->GetDofTypes();

        NodeBase* newNode;
        if (sameCoordinateNodes.original.node == nullptr)
        {
            // all nodeMap at this coordinate did not exist before
            // --> give random new ID
            int nodeId = rS.NodeCreateDOFs(dofTypes, coordinateNodePair.first);
            newNode = rS.NodeGetNodePtr(nodeId);
        }
        else
        {
            newNode = rS.NodePtrCreate(dofTypes, coordinateNodePair.first);
            int nodeId = nodeToId[sameCoordinateNodes.original.node];
            // replace the node ptr in all constraints, groups, and the original element
            rS.NodeExchangePtr(nodeId, sameCoordinateNodes.original.node, newNode, sameCoordinateNodes.original.elements);
        }

        // create one new node via the structure (to obtain a valid id)
        for (const auto& entry : sameCoordinateNodes.newEntries)
        {
            if (dofTypes != entry.node->GetDofTypes())
                throw MechanicsException(__PRETTY_FUNCTION__, "Dof types do not match... TODO improve msg.");
            entry.element->ExchangeNodePtr(entry.node, newNode);
            delete entry.node;
        }
    }
}
