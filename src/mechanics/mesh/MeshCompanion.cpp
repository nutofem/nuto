#include "mechanics/mesh/MeshCompanion.h"
#include "mechanics/structures/unstructured/Structure.h"

std::vector<NuTo::ElementBase*> NuTo::MeshCompanion::GetElementVector(Structure& rS, int rElementGroupId)
{
    std::vector<ElementBase*> elements;
    elements.reserve(rS.GroupGetNumMembers(rElementGroupId));

    for (int elementId : rS.GroupGetMemberIds(rElementGroupId))
        elements.push_back(rS.ElementGetElementPtr(elementId));

    return elements;
}

std::map<const NuTo::NodeBase*, int> NuTo::MeshCompanion::GetNodeToIdMap(NuTo::Structure& rS)
{
    std::map<const NodeBase*, int> nodeToId;
    int gAllNodes = rS.GroupGetNodesTotal();
    for (int nodeId : rS.GroupGetMemberIds(gAllNodes))
        nodeToId[rS.NodeGetNodePtr(nodeId)] = nodeId;
    rS.GroupDelete(gAllNodes);
    return nodeToId;
}