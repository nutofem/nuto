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