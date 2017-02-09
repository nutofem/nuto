#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/mesh/MeshCompanion.h"
#include "mechanics/structures/unstructured/Structure.h"

#include "mechanics/elements/ElementBase.h"
#include "mechanics/nodes/NodeBase.h"
#include "base/Timer.h"
#include "mechanics/interpolationtypes/InterpolationType.h"
#include "mechanics/interpolationtypes/InterpolationBase.h"

#include "math/SpacialContainer.h"

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





//! @brief store one newly created node and the corresponding element
//! @remark for replacing this node in the element
struct NodeElementPair
{
    NuTo::NodeBase* node;
    NuTo::ElementBase* element;
    bool existsBeforeConvert;
};

struct NodeElementPairCoordinate
{
    Eigen::VectorXd operator() (const NodeElementPair& rPair)
    {
        return rPair.node->Get(NuTo::Node::eDof::COORDINATES);
    }
};


std::vector<NodeElementPair> GetExistingNodes(NuTo::Structure& rS, std::vector<NuTo::ElementBase*> rElements)
{
    std::vector<NodeElementPair> existingNodes;
    for (NuTo::ElementBase* element : rElements)
    {
        // loop through nodes with coordinates
        for (int iNode = 0; iNode < element->GetNumNodes(NuTo::Node::eDof::COORDINATES); ++iNode)
        {
            NuTo::NodeBase* node = element->GetNode(iNode, NuTo::Node::eDof::COORDINATES);
            NodeElementPair pair;
            pair.node = node;
            pair.element = element;
            pair.existsBeforeConvert = true;
            existingNodes.push_back(pair);
        }
    }
    return existingNodes;
}

std::vector<NodeElementPair> CreateNewNodes(NuTo::Structure& rS, std::vector<NuTo::ElementBase*> rElements)
{
    std::vector<NodeElementPair> existingNodes;
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

            element->SetNode(iNode, node);

            NodeElementPair pair;
            pair.node = node;
            pair.element = element;
            pair.existsBeforeConvert = false;
            existingNodes.push_back(pair);
        }
    }
    return existingNodes;
}


std::set<NuTo::Node::eDof> GetDofTypeUnion(const std::vector<NodeElementPair>& pairs)
{
    std::set<NuTo::Node::eDof> allDofs;
    for (const auto& entry : pairs)
    {
        const auto& dofsOfNode = entry.node->GetDofTypes();
        std::set_union(allDofs.begin(), allDofs.end(),
                       dofsOfNode.begin(), dofsOfNode.end(),
                       std::inserter(allDofs, allDofs.end()));
    }
    return allDofs;
}



//! @brief returns the node that existed before calling ConvertToInterpolationType
//! @param rPairs vector of node-element-pairs
//! @return nodeptr of existing node or nullptr if rPairs to not contain a existing node
NuTo::NodeBase* GetExistingNode(const std::vector<NodeElementPair>& rPairs)
{
    for (const auto& pair : rPairs)
        if (pair.existsBeforeConvert)
            return pair.node;
    return nullptr;
}

std::vector<NuTo::ElementBase*> ExtractElementVector(const std::vector<NodeElementPair>& rPairs)
{
    std::vector<NuTo::ElementBase*> elements;
    elements.reserve(rPairs.size());
    for (const auto& pair : rPairs)
        elements.push_back(pair.element);
    return elements;
}

void NuTo::MeshCompanion::ElementConvertToInterpolationType(Structure& rS, int rGroupNumberElements, double rNodeDistanceMerge)
{
    Timer timer(__FUNCTION__, rS.GetShowTime(), rS.GetLogger());

    std::vector<ElementBase*> elements = GetElementVector(rS, rGroupNumberElements);

    std::vector<NodeElementPair> existingNodes = GetExistingNodes(rS, elements);
    std::vector<NodeElementPair> newNodes = CreateNewNodes(rS, elements);

    // concatenate...
    newNodes.insert(newNodes.end(), existingNodes.begin(), existingNodes.end());

    NuTo::SpacialContainer<NodeElementPair, NodeElementPairCoordinate> spacialContainer(newNodes);
    auto duplicates = spacialContainer.GetAllDuplicateValues(rNodeDistanceMerge);

    std::cout << "num duplicates " << duplicates.size() << std::endl;

    std::map<const NodeBase*, int> nodeToId = GetNodeToIdMap(rS);
    for (auto& pairsAtSameCoordinate : duplicates)
    {
        Eigen::VectorXd coordinate = NodeElementPairCoordinate()(pairsAtSameCoordinate[0]);
        auto dofTypes = GetDofTypeUnion(pairsAtSameCoordinate);

        NuTo::NodeBase* newNode;

        NuTo::NodeBase* existingNode = GetExistingNode(pairsAtSameCoordinate);
        if (existingNode)
        {
            int id = nodeToId[existingNode];
            newNode = rS.NodePtrCreate(dofTypes, coordinate);

            // replace the node ptr in all constraints, groups, and the original element
            rS.NodeExchangePtr(id, existingNode, newNode, {}); // deletes existingNode
        }
        else
        {
            int id = rS.NodeCreate(coordinate, dofTypes);
            newNode = rS.NodeGetNodePtr(id);
        }

        for (auto& pair : pairsAtSameCoordinate)
        {
            pair.element->ExchangeNodePtr(pair.node, newNode);
            if (not pair.existsBeforeConvert)
                delete pair.node;
        }
    }
}

