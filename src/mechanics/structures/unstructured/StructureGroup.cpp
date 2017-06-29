#include "mechanics/structures/unstructured/Structure.h"

#include "base/Timer.h"
#include "mechanics/elements/ElementBase.h"
#include "mechanics/groups/Group.h"
#include "mechanics/groups/GroupEnum.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"

#include <eigen3/Eigen/Dense>

void NuTo::Structure::GroupAddElementFromType(int rIdentGroup, int rInterpolationType)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());
    boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rIdentGroup);
    if (itGroup==mGroupMap.end())
        throw Exception(__PRETTY_FUNCTION__, "Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=eGroupId::Elements)
        throw Exception(__PRETTY_FUNCTION__, "Group is not an element group.");

    boost::ptr_map<int,InterpolationType>::iterator itInterpolationType = mInterpolationTypeMap.find(rInterpolationType);
    if (itInterpolationType==mInterpolationTypeMap.end())
        throw Exception(__PRETTY_FUNCTION__, "InterpolationType with the given identifier does not exist.");

    InterpolationType* interpolationType = itInterpolationType->second;

    for (boost::ptr_map<int,ElementBase>::iterator elementIter = this->mElementMap.begin(); elementIter != this->mElementMap.end(); elementIter++)
    {
        ElementBase* element = elementIter->second;

        if (&element->GetInterpolationType() == interpolationType)
            itGroup->second->AddMember(elementIter->first, elementIter->second);
    }
}


void NuTo::Structure::GroupAddElement(int rIdentGroup, int rIdElement)
{
    boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rIdentGroup);
    if (itGroup==mGroupMap.end())
        throw Exception(__PRETTY_FUNCTION__, "Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=eGroupId::Elements)
        throw Exception(__PRETTY_FUNCTION__, "An element can be added only to an element group.");

    itGroup->second->AddMember(rIdElement, ElementGetElementPtr(rIdElement));
}

void NuTo::Structure::GroupAddElementsTotal(int rIdentGroup)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rIdentGroup);

    if (itGroup==mGroupMap.end())
        throw Exception(__PRETTY_FUNCTION__, "Group with the given identifier does not exist.");

    if (itGroup->second->GetType()!=eGroupId::Elements)
        throw Exception(__PRETTY_FUNCTION__, "An element can be added only to an element group.");

    for (auto const& iPair : mElementMap)
        itGroup->second->AddMember(iPair.first, iPair.second);
}

int NuTo::Structure::GroupGetElementsTotal()
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    int groupId = GroupCreate(eGroupId::Elements);
    auto& elementGroup = *(mGroupMap.at(groupId).AsGroupElement());
    for (auto const& iPair : mElementMap)
        elementGroup.AddMember(iPair.first, iPair.second);
    return groupId;
}

int NuTo::Structure::GroupGetNodesTotal()
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    int groupId = GroupCreate(eGroupId::Nodes);
    auto& elementGroup = *(mGroupMap.at(groupId).AsGroupNode());
    for (auto const& iPair : mNodeMap)
        elementGroup.AddMember(iPair.first, iPair.second);
    return groupId;
}

void NuTo::Structure::GroupAddNodeFromElementGroupCoordinateRange(int rIdentNodeGroup, int rSearchIdentElementGroup, int rDirection, double rMin, double rMax)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    boost::ptr_map<int, GroupBase>::iterator itGroup = mGroupMap.find(rIdentNodeGroup);
    if (itGroup == mGroupMap.end())
        throw Exception(__PRETTY_FUNCTION__, "Group with the given identifier does not exist.");
    if (itGroup->second->GetType() != eGroupId::Nodes)
        throw Exception(__PRETTY_FUNCTION__, "A node can be added only to a node group.");

    if (rDirection < 0 || rDirection > mDimension)
        throw Exception(__PRETTY_FUNCTION__, "The direction is either 0(x),1(Y) or 2(Z) and has to be smaller than the dimension of the structure.");

    auto elementsInGroup = GroupGetMemberIds(rSearchIdentElementGroup);

    std::set<int> nodesInGroup;
    for (auto elementId : elementsInGroup)
        for (int nodeId : ElementGetNodes(elementId))
            nodesInGroup.insert(nodeId);


    for (int iNodeId : nodesInGroup)
    {
        auto nodePtr = NodeGetNodePtr(iNodeId);
        double coordinate = nodePtr->Get(Node::eDof::COORDINATES)[rDirection];

        if (coordinate >= rMin and coordinate <= rMax)
            itGroup->second->AddMember(iNodeId, nodePtr);
    }
}

