#include "mechanics/mesh/MeshCompanion.h"
#include "mechanics/structures/unstructured/Structure.h"

std::vector<NuTo::ElementBase*> NuTo::MeshCompanion::GetElementVector(Structure& s, int elementGroupId)
{
    std::vector<ElementBase*> elements;
    elements.reserve(s.GroupGetNumMembers(elementGroupId));

    for (int elementId : s.GroupGetMemberIds(elementGroupId))
        elements.push_back(s.ElementGetElementPtr(elementId));

    return elements;
}

std::map<const NuTo::NodeBase*, int> NuTo::MeshCompanion::GetNodeToIdMap(NuTo::Structure& s)
{
    std::map<const NodeBase*, int> nodeToId;
    int gAllNodes = s.GroupGetNodesTotal();
    for (int nodeId : s.GroupGetMemberIds(gAllNodes))
        nodeToId[s.NodeGetNodePtr(nodeId)] = nodeId;
    s.GroupDelete(gAllNodes);
    return nodeToId;
}
