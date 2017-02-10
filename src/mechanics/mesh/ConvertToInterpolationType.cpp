#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/mesh/MeshCompanion.h"
#include "mechanics/structures/unstructured/Structure.h"

#include "mechanics/elements/ElementBase.h"
#include "mechanics/nodes/NodeBase.h"
#include "base/Timer.h"
#include "mechanics/interpolationtypes/InterpolationType.h"
#include "mechanics/interpolationtypes/InterpolationBase.h"

#include "math/SpacialContainer.h"
#include <limits>


//! @brief calculate a reasonable node merge distance
double GetMergeDistance(NuTo::Structure& rS, int rGroupNumberElements)
{
    // calculate and store the 'size' of each element ...
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
    // length ~ area**(1/2) ~ volume**(1/3)
    double lengthMin = std::pow(smallestElementSize, 1. / rS.GetDimension());
    return std::max(lengthMin / 1000., std::numeric_limits<double>::epsilon() * 10);
}

void NuTo::MeshCompanion::ElementTotalConvertToInterpolationType(Structure& rS)
{
    int groupNumber = rS.GroupGetElementsTotal();
    ElementConvertToInterpolationType(rS, groupNumber);
    rS.GroupDelete(groupNumber);
}

void NuTo::MeshCompanion::ElementConvertToInterpolationType(Structure& rS, int rGroupNumberElements)
{
    ElementConvertToInterpolationType(rS, rGroupNumberElements, GetMergeDistance(rS, rGroupNumberElements));
}


void NuTo::MeshCompanion::ElementTotalConvertToInterpolationType(Structure& rS, double rNodeDistanceMerge)
{
    int groupNumber = rS.GroupGetElementsTotal();
    ElementConvertToInterpolationType(rS, groupNumber, rNodeDistanceMerge);
    rS.GroupDelete(groupNumber);
}

static constexpr int NOT_SET = -1337;

//! @brief node information required to build
struct TmpNode
{
    Eigen::VectorXd coordinate;
    std::set<NuTo::Node::eDof> dofs;
    NuTo::ElementBase* element;
    int elementNodeId;
    int originalId = NOT_SET;
};

//! @brief struct to extract the coodinates from TmpNode (for NuTo::SpacialContainer)
struct TmpNodeCoordinate
{
    Eigen::VectorXd operator() (const TmpNode& rTmpNode)
    {
        return rTmpNode.coordinate;
    }
};

//! @brief extract TmpNodes from existing nodes
std::vector<TmpNode> GetExistingTmpNodes(NuTo::Structure& rS, std::vector<NuTo::ElementBase*> rElements, std::map<const NuTo::NodeBase*, int> rNodeToId)
{
    std::vector<TmpNode> existingNodes;
    for (NuTo::ElementBase* element : rElements)
    {
        // loop through nodes with coordinates
        for (int iNode = 0; iNode < element->GetNumNodes(NuTo::Node::eDof::COORDINATES); ++iNode)
        {
            NuTo::NodeBase* node = element->GetNode(iNode, NuTo::Node::eDof::COORDINATES);
            TmpNode pair;
            pair.coordinate = node->Get(NuTo::Node::eDof::COORDINATES);
            pair.dofs = {};
            pair.element = element;
            pair.elementNodeId = element->GetInterpolationType().Get(NuTo::Node::eDof::COORDINATES).GetNodeIndex(iNode);
            pair.originalId = rNodeToId[node];

            existingNodes.push_back(pair);
        }
    }
    return existingNodes;
}

//! @brief create TmpNodes from existing nodes
std::vector<TmpNode> GetNewTmpNodes(NuTo::Structure& rS, std::vector<NuTo::ElementBase*> rElements)
{
    std::vector<TmpNode> existingNodes;
    for (NuTo::ElementBase* element : rElements)
    {
        const NuTo::InterpolationType& interpolationType = element->GetInterpolationType();
        element->ResizeNodes(interpolationType.GetNumNodes());
        for (int iNode = 0; iNode < element->GetNumNodes(); ++iNode)
        {
            Eigen::VectorXd naturalNodeCoordinates = interpolationType.GetNaturalNodeCoordinates(iNode);
            Eigen::VectorXd globalNodeCoordinates = element->InterpolateDofGlobal(naturalNodeCoordinates, NuTo::Node::eDof::COORDINATES);

            TmpNode pair;
            pair.coordinate = globalNodeCoordinates;
            pair.dofs = interpolationType.GetNodeDofs(iNode);
            pair.element = element;
            pair.elementNodeId = iNode;
            existingNodes.push_back(pair);
        }
    }
    return existingNodes;
}

//! @brief unites the dof types of all tmp nodes at the same coordinate
std::set<NuTo::Node::eDof> GetDofTypeUnion(const std::vector<TmpNode>& rSameCoordinateNodes)
{
    std::set<NuTo::Node::eDof> allDofs;
    for (const auto& entry : rSameCoordinateNodes)
    {
        std::set_union(allDofs.begin(), allDofs.end(),
                       entry.dofs.begin(), entry.dofs.end(),
                       std::inserter(allDofs, allDofs.end()));
    }
    return allDofs;
}

//! @brief finds the old node id (before calling ElementConvertToInterpolationType)
//! @return node id, or NOT_SET if the node is new
int FindPreviousNodeId(const std::vector<TmpNode>& rSameCoordinateNodes)
{
    for (const auto& node : rSameCoordinateNodes)
        if (node.originalId != NOT_SET)
            return node.originalId;
    return NOT_SET;
}

void NuTo::MeshCompanion::ElementConvertToInterpolationType(Structure& rS, int rGroupNumberElements, double rNodeDistanceMerge)
{
    Timer timer(__FUNCTION__, rS.GetShowTime(), rS.GetLogger());

    std::vector<ElementBase*> elements = GetElementVector(rS, rGroupNumberElements);

    std::vector<TmpNode> existingNodes = GetExistingTmpNodes(rS, elements, GetNodeToIdMap(rS));
    std::vector<TmpNode> newNodes = GetNewTmpNodes(rS, elements);

    // concatenate...
    newNodes.insert(newNodes.end(), existingNodes.begin(), existingNodes.end());

    NuTo::SpacialContainer<TmpNode, TmpNodeCoordinate> spacialContainer(newNodes);
    auto duplicates = spacialContainer.GetAllDuplicateValues(rNodeDistanceMerge);

    for (auto& pairsAtSameCoordinate : duplicates)
    {
        Eigen::VectorXd coordinate = pairsAtSameCoordinate[0].coordinate;
        auto dofTypes = GetDofTypeUnion(pairsAtSameCoordinate);

        int oldNodeId = FindPreviousNodeId(pairsAtSameCoordinate);
        if (oldNodeId == NOT_SET)
        {
            NuTo::NodeBase* newNode = rS.NodeGetNodePtr(rS.NodeCreate(coordinate, dofTypes));
            for (auto& pair : pairsAtSameCoordinate)
                pair.element->SetNode(pair.elementNodeId, newNode);
        }
        else
        {
            NuTo::NodeBase* newNode = rS.NodePtrCreate(dofTypes, coordinate);
            std::vector<ElementBase*> elementsToChange;
            for (auto& pair : pairsAtSameCoordinate)
                elementsToChange.push_back(pair.element);

            // replace the node ptr in all constraints, groups, and the original element
            rS.NodeExchangePtr(oldNodeId, rS.NodeGetNodePtr(oldNodeId), newNode, elementsToChange);
        }
    }
}

