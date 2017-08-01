#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/mesh/MeshCompanion.h"
#include "mechanics/structures/unstructured/Structure.h"

#include "mechanics/elements/ElementBase.h"
#include "mechanics/nodes/NodeBase.h"
#include "base/Timer.h"
#include "mechanics/interpolationtypes/InterpolationType.h"
#include "mechanics/interpolationtypes/InterpolationBase.h"

#include "math/SpatialContainer.h"
#include <limits>


//! @brief calculate a reasonable node merge distance
double GetMergeDistance(NuTo::Structure& s, int elementGroup)
{
    // calculate and store the 'size' of each element ...
    // = volume in 3D
    // = area in 2D
    // = length in 1D
    double smallestElementSize = 1e42;
    for (int elementId : s.GroupGetMemberIds(elementGroup))
    {
        NuTo::ElementBase* element = s.ElementGetElementPtr(elementId);
        Eigen::VectorXd sizeForEachIntegrationPoint = element->GetIntegrationPointVolume();
        smallestElementSize = std::min(smallestElementSize, sizeForEachIntegrationPoint.sum());
    }
    // length ~ area**(1/2) ~ volume**(1/3)
    double lengthMin = std::pow(smallestElementSize, 1. / s.GetDimension());
    return std::max(lengthMin / 1000., std::numeric_limits<double>::epsilon() * 10);
}

void NuTo::MeshCompanion::ElementTotalConvertToInterpolationType(Structure& s)
{
    int groupNumber = s.GroupGetElementsTotal();
    ElementConvertToInterpolationType(s, groupNumber);
    s.GroupDelete(groupNumber);
}

void NuTo::MeshCompanion::ElementConvertToInterpolationType(Structure& s, int elementGroup)
{
    ElementConvertToInterpolationType(s, elementGroup, GetMergeDistance(s, elementGroup));
}


void NuTo::MeshCompanion::ElementTotalConvertToInterpolationType(Structure& s, double nodeMergeDistance)
{
    int groupNumber = s.GroupGetElementsTotal();
    ElementConvertToInterpolationType(s, groupNumber, nodeMergeDistance);
    s.GroupDelete(groupNumber);
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

//! @brief struct to extract the coodinates from TmpNode (for NuTo::SpatialContainer)
struct TmpNodeCoordinate
{
    Eigen::VectorXd operator()(const TmpNode& tmpNode) const
    {
        return tmpNode.coordinate;
    }
};

//! @brief extract TmpNodes from existing nodes
std::vector<TmpNode> GetExistingTmpNodes(NuTo::Structure&, std::vector<NuTo::ElementBase*> elements,
                                         std::map<const NuTo::NodeBase*, int> nodeToId)
{
    std::vector<TmpNode> existingNodes;
    for (NuTo::ElementBase* element : elements)
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
            pair.originalId = nodeToId[node];

            existingNodes.push_back(pair);
        }
    }
    return existingNodes;
}

//! @brief create TmpNodes from existing nodes
std::vector<TmpNode> GetNewTmpNodes(NuTo::Structure&, std::vector<NuTo::ElementBase*> elements)
{
    std::vector<TmpNode> existingNodes;
    for (NuTo::ElementBase* element : elements)
    {
        const NuTo::InterpolationType& interpolationType = element->GetInterpolationType();
        element->ResizeNodes(interpolationType.GetNumNodes());
        for (int iNode = 0; iNode < element->GetNumNodes(); ++iNode)
        {
            Eigen::VectorXd naturalNodeCoordinates = interpolationType.GetNaturalNodeCoordinates(iNode);
            Eigen::VectorXd globalNodeCoordinates =
                    element->InterpolateDofGlobal(naturalNodeCoordinates, NuTo::Node::eDof::COORDINATES);

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
std::set<NuTo::Node::eDof> GetDofTypeUnion(const std::vector<TmpNode>& sameCoordinateNodes)
{
    std::set<NuTo::Node::eDof> allDofs;
    for (const auto& entry : sameCoordinateNodes)
    {
        std::set_union(allDofs.begin(), allDofs.end(), entry.dofs.begin(), entry.dofs.end(),
                       std::inserter(allDofs, allDofs.end()));
    }
    return allDofs;
}

//! @brief finds the old node id (before calling ElementConvertToInterpolationType)
//! @return node id, or NOT_SET if the node is new
int FindPreviousNodeId(const std::vector<TmpNode>& sameCoordinateNodes)
{
    for (const auto& node : sameCoordinateNodes)
        if (node.originalId != NOT_SET)
            return node.originalId;
    return NOT_SET;
}

void NuTo::MeshCompanion::ElementConvertToInterpolationType(Structure& s, int elementGroup,
                                                            double nodeMergeDistance)
{
    Timer timer(__FUNCTION__, s.GetShowTime(), s.GetLogger());

    std::vector<ElementBase*> elements = GetElementVector(s, elementGroup);

    std::vector<TmpNode> existingNodes = GetExistingTmpNodes(s, elements, GetNodeToIdMap(s));
    std::vector<TmpNode> newNodes = GetNewTmpNodes(s, elements);

    // concatenate...
    newNodes.insert(newNodes.end(), existingNodes.begin(), existingNodes.end());

    NuTo::SpatialContainer<TmpNode, TmpNodeCoordinate> spacialContainer(newNodes);
    auto duplicates = spacialContainer.GetAllDuplicateValues(nodeMergeDistance);

    for (auto& pairsAtSameCoordinate : duplicates)
    {
        Eigen::VectorXd coordinate = pairsAtSameCoordinate[0].coordinate;
        auto dofTypes = GetDofTypeUnion(pairsAtSameCoordinate);

        int oldNodeId = FindPreviousNodeId(pairsAtSameCoordinate);
        if (oldNodeId == NOT_SET)
        {
            NuTo::NodeBase* newNode = s.NodeGetNodePtr(s.NodeCreate(coordinate, dofTypes));
            for (auto& pair : pairsAtSameCoordinate)
                pair.element->SetNode(pair.elementNodeId, newNode);
        }
        else
        {
            NuTo::NodeBase* newNode = s.NodePtrCreate(dofTypes, coordinate);
            std::vector<ElementBase*> elementsToChange;
            for (auto& pair : pairsAtSameCoordinate)
                elementsToChange.push_back(pair.element);

            // replace the node ptr in all constraints, groups, and the original element
            s.NodeExchangePtr(oldNodeId, s.NodeGetNodePtr(oldNodeId), newNode, elementsToChange);
        }
    }
}
