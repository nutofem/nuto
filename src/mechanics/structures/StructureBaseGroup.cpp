#include <set>
#include "base/Timer.h"

#include "mechanics/structures/StructureBase.h"
#include "mechanics/groups/Group.h"
#include "mechanics/groups/GroupEnum.h"
#include "mechanics/elements/ElementBase.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/DirectionEnum.h"

void NuTo::StructureBase::GroupInfo(int rVerboseLevel) const
{
    NuTo::Timer timer(__FUNCTION__, GetShowTime(), GetLogger());
    std::cout << "number of groups  : " << mGroupMap.size() << std::endl;
    if (rVerboseLevel > 2)
    {
        for (boost::ptr_map<int, GroupBase>::const_iterator it = mGroupMap.begin(); it != mGroupMap.end(); it++)
        {
            std::cout << "  Group " << it->first << std::endl;
            it->second->Info(rVerboseLevel, this);
            std::cout << std::endl;
        }
    }
}

int NuTo::StructureBase::GroupGetId(GroupBase* rGroup) const
{
    for (boost::ptr_map<int, GroupBase>::const_iterator it = mGroupMap.begin(); it != mGroupMap.end(); it++)
    {
        if (it->second == rGroup)
            return it->first;
    }
    throw MechanicsException("[NuTo::Structure::GroupGetId] Group does not exist.");
}

int NuTo::StructureBase::GroupCreate(const std::string& rType)
{
    NuTo::Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    //find unused integer id
    int groupNumber = GetUnusedId(mGroupMap);

    // transform string to uppercase
    std::string GroupTypeString;
    std::transform(rType.begin(), rType.end(), std::back_inserter(GroupTypeString), (int (*)(int)) toupper);

    if(    GroupTypeString==std::string("ELEMENTS"))
    {
        GroupCreate(groupNumber,NuTo::eGroupId::Elements);
    }
    else
    {
        if (GroupTypeString==std::string("NODES"))
        {
            GroupCreate(groupNumber,NuTo::eGroupId::Nodes);
        }
        else
        {
            throw MechanicsException("[NuTo::StructureBase::GroupCreate] Group type not implemented.");
        }
    }
    return groupNumber;
}

int NuTo::StructureBase::GroupCreate(NuTo::eGroupId rEnumType)
{
    NuTo::Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    int groupNumber = GetUnusedId(mGroupMap);
    GroupCreate(groupNumber, rEnumType);
    return groupNumber;
}

void NuTo::StructureBase::GroupCreate(int id, NuTo::eGroupId rEnumType)
{
    boost::ptr_map<int, GroupBase>::iterator it = mGroupMap.find(id);
    if (it != mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupCreate] Group id already exists.");

    switch (rEnumType)
    {
    case NuTo::eGroupId::Elements:
        mGroupMap.insert(id, new NuTo::Group<ElementBase>);
        break;
    case NuTo::eGroupId::Nodes:
        mGroupMap.insert(id, new NuTo::Group<NodeBase>);
        break;
    default:
        throw MechanicsException("[NuTo::StructureBase::GroupCreate] Group type not implemented.");
    }
}

void NuTo::StructureBase::GroupDelete(int rIdentGroup)
{
    NuTo::Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    // find group in map
    boost::ptr_map<int, GroupBase>::iterator itGroup = mGroupMap.find(rIdentGroup);
    if (itGroup == this->mGroupMap.end())
    {
        throw MechanicsException("[NuTo::StructureBase::GroupDelete] Group with the given identifier does not exist.");
    }
    else
    {
        // delete load from map
        this->mGroupMap.erase(itGroup);
    }
}

void NuTo::StructureBase::GroupAddNode(int rIdentGroup, int rIdNode)
{
    boost::ptr_map<int, GroupBase>::iterator itGroup = mGroupMap.find(rIdentGroup);
    if (itGroup == mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupAddNode] Group with the given identifier does not exist.");
    if (itGroup->second->GetType() != eGroupId::Nodes)
        throw MechanicsException("[NuTo::StructureBase::GroupAddNode] A node can be added only to a node group.");

    itGroup->second->AddMember(rIdNode, NodeGetNodePtr(rIdNode));
}

void NuTo::StructureBase::GroupAddElement(int rIdentGroup, int rIdElement)
{
    boost::ptr_map<int, GroupBase>::iterator itGroup = mGroupMap.find(rIdentGroup);
    if (itGroup == mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupAddElement] Group with the given identifier does not exist.");
    if (itGroup->second->GetType() != eGroupId::Elements)
        throw MechanicsException("[NuTo::StructureBase::GroupAddElement] An element can be added only to an element group.");

    itGroup->second->AddMember(rIdElement, ElementGetElementPtr(rIdElement));
}

void NuTo::StructureBase::GroupAddNodeCoordinateRange(int rIdentGroup, int rDirection, double rMin, double rMax)
{
    NuTo::Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    boost::ptr_map<int, GroupBase>::iterator itGroup = mGroupMap.find(rIdentGroup);
    if (itGroup == mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupAddNodeCoordinateRange] Group with the given identifier does not exist.");
    if (itGroup->second->GetType() != eGroupId::Nodes)
        throw MechanicsException("[NuTo::StructureBase::GroupAddNodeCoordinateRange] A node can be added only to a node group.");

    if (rDirection < 0 || rDirection > mDimension)
        throw MechanicsException("[NuTo::StructureBase::GroupAddNodeCoordinateRange] The direction is either 0(x),1(Y) or 2(Z) and has to be smaller than the dimension of the structure.");

    std::vector<std::pair<int, NodeBase*> > nodeVector;
    this->GetNodesTotal(nodeVector);

    for (auto& node : nodeVector)
    {
        NodeBase* nodePtr(node.second);
        if (nodePtr->GetNum(Node::eDof::COORDINATES) < 1)
            continue;
        double coordinate = nodePtr->Get(Node::eDof::COORDINATES)[rDirection];

        if (coordinate >= rMin && coordinate <= rMax)
            itGroup->second->AddMember(node.first, nodePtr);
    }
}

NuTo::Group<NuTo::NodeBase>& NuTo::StructureBase::GroupGetNodeCoordinateRange(eDirection direction, double min, double max)
{
    int groupId = GroupCreate(eGroupId::Nodes);
    GroupAddNodeCoordinateRange(groupId, direction, min, max);
    return *GroupGetGroupPtr(groupId)->AsGroupNode();  
}

NuTo::Group<NuTo::NodeBase>& NuTo::StructureBase::GroupGetNodesAtCoordinate(eDirection direction, double value, double tolerance)
{
    return GroupGetNodeCoordinateRange(direction, value - tolerance, value + tolerance);
}

void NuTo::StructureBase::GroupAddNodeCoordinateRange(int rIdentGroup, NuTo::eDirection rDirection, double rMin, double rMax)
{
    GroupAddNodeCoordinateRange(rIdentGroup, ToComponentIndex(rDirection), rMin, rMax);
}

void NuTo::StructureBase::GroupAddNodeFunction(int rIdentNewGroup, int rIdentOldGroup,  std::function<bool(NuTo::NodeBase *)> rFunction)
{
    NuTo::Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    boost::ptr_map<int, GroupBase>::iterator itGroupOld = mGroupMap.find(rIdentOldGroup);
    if (itGroupOld == mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupAddNodeCoordinateRange] Old group with the given identifier does not exist.");
    if (itGroupOld->second->GetType() != eGroupId::Nodes)
        throw MechanicsException("[NuTo::StructureBase::GroupAddNodeCoordinateRange] A node can be added only to a node group.");

    boost::ptr_map<int, GroupBase>::iterator itGroupNew = mGroupMap.find(rIdentNewGroup);
    if (itGroupNew == mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupAddNodeCoordinateRange] New group with the given identifier does not exist.");
    if (itGroupNew->second->GetType() != eGroupId::Nodes)
        throw MechanicsException("[NuTo::StructureBase::GroupAddNodeCoordinateRange] A node can be added only to a node group.");


    std::vector<int> members;
    this->NodeGroupGetMembers(rIdentOldGroup,  members);

    for (int memberId : members)
    {
        NodeBase* nodePtr = this->NodeGetNodePtr(memberId);
        if (nodePtr->GetNum(Node::eDof::COORDINATES) < 1)
            continue;

        if (rFunction(nodePtr))
            itGroupNew->second->AddMember(memberId, nodePtr);
    }
}

void NuTo::StructureBase::GroupAddNodeFunction(int rIdentGroup, std::function<bool(NuTo::NodeBase *)> rFunction)
{
    NuTo::Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    boost::ptr_map<int, GroupBase>::iterator itGroup = mGroupMap.find(rIdentGroup);
    if (itGroup == mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupAddNodeFunction] Group with the given identifier does not exist.");
    if (itGroup->second->GetType() != eGroupId::Nodes)
        throw MechanicsException("[NuTo::StructureBase::GroupAddNodeFunction] A node can be added only to a node group.");

    std::vector<std::pair<int, NodeBase*> > nodeVector;
    this->GetNodesTotal(nodeVector);

    for (auto& node : nodeVector)
    {
        NodeBase* nodePtr(node.second);
        if (rFunction(nodePtr))
            itGroup->second->AddMember(node.first, nodePtr);
    }
}

void NuTo::StructureBase::GroupAddNodeRadiusRange(int rIdentGroup, Eigen::VectorXd rCenter, double rMin, double rMax)
{
    NuTo::Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    boost::ptr_map<int, GroupBase>::iterator itGroup = mGroupMap.find(rIdentGroup);
    if (itGroup == mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupAddNodeRadiusRange] Group with the given identifier does not exist.");
    if (itGroup->second->GetType() != eGroupId::Nodes)
        throw MechanicsException("[NuTo::StructureBase::GroupAddNodeRadiusRange] A node can be added only to a node group.");

    if (rCenter.rows() != mDimension || rCenter.cols() != 1)
        throw MechanicsException("[NuTo::StructureBase::GroupAddNodeRadiusRange] The center point must have the same number of coordinates as the dimension of the structure.");

    if (rMin > rMax)
        throw MechanicsException("[NuTo::StructureBase::GroupAddNodeRadiusRange] The minimum radius must not be larger than the maximum radius.");

    std::vector<std::pair<int, NodeBase*> > nodeVector;
    this->GetNodesTotal(nodeVector);
    double rMin2 = rMin * rMin;
    double rMax2 = rMax * rMax;

    for (auto& node : nodeVector)
    {
        NodeBase* nodePtr(node.second);
        if (nodePtr->GetNum(Node::eDof::COORDINATES) < 1)
            continue;
        Eigen::VectorXd dCoordinates = nodePtr->Get(Node::eDof::COORDINATES) - rCenter;
        double r2 = dCoordinates.dot(dCoordinates);

        if (r2 >= rMin2 && r2 <= rMax2)
            itGroup->second->AddMember(node.first, nodePtr);
    }
}

void NuTo::StructureBase::GroupAddNodeCylinderRadiusRange(int rIdentGroup, Eigen::VectorXd rCenter, Eigen::VectorXd rDirection, double rMin, double rMax)
{
    NuTo::Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    boost::ptr_map<int, GroupBase>::iterator itGroup = mGroupMap.find(rIdentGroup);
    if (itGroup == mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupAddNodeCylinderRadiusRange] Group with the given identifier does not exist.");
    if (itGroup->second->GetType() != eGroupId::Nodes)
        throw MechanicsException("[NuTo::StructureBase::GroupAddNodeCylinderRadiusRange] A node can be added only to a node group.");

    if (rCenter.rows() != mDimension || rCenter.cols() != 1)
        throw MechanicsException("[NuTo::StructureBase::GroupAddNodeCylinderRadiusRange] The center point must have the same number of coordinates as the dimension of the structure.");

    if (rDirection.rows() != mDimension || rDirection.cols() != 1)
        throw MechanicsException("[NuTo::StructureBase::GroupAddNodeCylinderRadiusRange] The direction point must have the same number of coordinates as the dimension of the structure.");

    if (rMin > rMax)
        throw MechanicsException("[NuTo::StructureBase::GroupAddNodeCylinderRadiusRange] The minimum radius must not be larger than the maximum radius.");

    switch (mDimension)
    {
    case 2:
    {
        std::vector<std::pair<int, NodeBase*> > nodeVector;
        this->GetNodesTotal(nodeVector);
        Eigen::Vector2d coordinates;
        Eigen::Vector2d vecDelta;
        double rMin2 = rMin * rMin;
        double rMax2 = rMax * rMax;

        for (auto& node : nodeVector)
        {
            NodeBase* nodePtr(node.second);
            if (nodePtr->GetNum(Node::eDof::COORDINATES) != 2)
                continue;
            double r2(0.);
            coordinates = nodePtr->Get(Node::eDof::COORDINATES);
            vecDelta = coordinates - rCenter;

            r2 = (vecDelta(0) * vecDelta(0) + vecDelta(1) * vecDelta(1));

            if (r2 >= rMin2 && r2 <= rMax2)
            {
                itGroup->second->AddMember(node.first, nodePtr);
            }
        }
        break;
    }
    case 3:
    {
        std::vector<std::pair<int, NodeBase*> > nodeVector;
        this->GetNodesTotal(nodeVector);
        Eigen::Vector3d coordinates;
        Eigen::Vector3d vecPtrCenter;
        Eigen::Vector3d vecPtrProjection;
        Eigen::Vector3d vecDelta;
        double rMin2 = rMin * rMin;
        double rMax2 = rMax * rMax;

        //normalize Diretion Vector
        rDirection *= 1. / rDirection.norm();

        for (auto& node : nodeVector)
        {
            NodeBase* nodePtr(node.second);
            if (nodePtr->GetNum(Node::eDof::COORDINATES) != 3)
                continue;
            double r2(0.);
            coordinates = nodePtr->Get(Node::eDof::COORDINATES);
            vecPtrCenter = coordinates - rCenter;

            // get projection onto axis
            double s = rDirection.transpose() * vecPtrCenter;
            vecPtrProjection = rCenter + rDirection * s;
            vecDelta = coordinates - vecPtrProjection;

            r2 = (vecDelta(0) * vecDelta(0) + vecDelta(1) * vecDelta(1) + vecDelta(2) * vecDelta(2));

            if (r2 >= rMin2 && r2 <= rMax2)
            {
                itGroup->second->AddMember(node.first, nodePtr);
            }
        }
        break;
    }
    default:
        throw MechanicsException("[NuTo::StructureBase::GroupAddNodeCylinderRadiusRange] Only implemented for 2D and 3D.");

    }
}

void NuTo::StructureBase::GroupAddElementsFromNodes(int rElementGroupId, int rNodeGroupId, bool rHaveAllNodes)
{
    NuTo::Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    //get element group
    Group<ElementBase> *elementGroup = this->GroupGetGroupPtr(rElementGroupId)->AsGroupElement();

    //get node group
    const Group<NodeBase> *nodeGroup = this->GroupGetGroupPtr(rNodeGroupId)->AsGroupNode();

    //since the search is done via the id's, the element nodes are ptr, so make another set with the node ptrs
    std::set<const NodeBase*> nodePtrSet;
    for (auto node : *nodeGroup)
    {
        nodePtrSet.insert(node.second);
    }

    std::vector<std::pair<int, ElementBase*> > elementVector;
    this->GetElementsTotal(elementVector);
    std::vector<const NodeBase*> elementNodes;
    for (auto& element : elementVector)
    {
        if (!elementGroup->Contain(element.first))
        {
            bool addElement;
            if (rHaveAllNodes == true)
            {
                addElement = true;
                int countNode;
                for (countNode = 0; (countNode < element.second->GetNumNodes()) && (addElement == true); countNode++)
                {
                    if (nodePtrSet.find(element.second->GetNode(countNode)) == nodePtrSet.end())
                    {
                        addElement = false;
                    }
                }
            }
            else
            {
                addElement = false;
                int countNode;
                for (countNode = 0; (countNode < element.second->GetNumNodes()) && (addElement == false); countNode++)
                {
                    if (nodePtrSet.find(element.second->GetNode(countNode)) != nodePtrSet.end())
                    {
                        addElement = true;
                    }
                }
            }
            if (addElement)
            {
                // add the element;
                elementGroup->AddMember(element.first, element.second);
            }
        }
    }
}

void NuTo::StructureBase::GroupAddNodesFromElements(int rNodeGroupId, int rElementGroupId)
{
    NuTo::Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    //get element group
    Group<ElementBase> *elementGroup = this->GroupGetGroupPtr(rElementGroupId)->AsGroupElement();

    //get node group
    Group<NodeBase> *nodeGroup = this->GroupGetGroupPtr(rNodeGroupId)->AsGroupNode();


    // create a map for fast node index search
    std::vector<std::pair<int,const NodeBase*> > nodesAndIds;
    GetNodesTotal(nodesAndIds);
    std::map<const NodeBase*, int> nodeToId;
    for (auto pair : nodesAndIds)
    {
        nodeToId[pair.second] = pair.first;
    }

    for (Group<ElementBase>::const_iterator itElement = elementGroup->begin(); itElement != elementGroup->end(); itElement++)
    {
        try
        {
            ElementBase* element = itElement->second;
            for (int iNode = 0; iNode < element->GetNumNodes(); ++iNode)
            {
                NuTo::NodeBase* nodePtr = element->GetNode(iNode);
                int nodeId = nodeToId[nodePtr];

                if (not nodeGroup->Contain(nodeId))
                    nodeGroup->AddMember(nodeId, nodePtr);
            }
        }
        catch (NuTo::MechanicsException &e)
        {
            e.AddMessage("[NuTo::StructureBase::GroupAddNodesFromElements] Error for element " + std::to_string(ElementGetId(itElement->second)) + ".");
            throw;
        }
    }
}


int NuTo::StructureBase::GroupCreateNodeGroupFromElements(int rElementGroupId)
{
    int n = GroupCreate(NuTo::eGroupId::Nodes);
    GroupAddNodesFromElements(n, rElementGroupId);
    return n;
}


int NuTo::StructureBase::GroupUnion(int rIdentGroup1, int rIdentGroup2)
{
    NuTo::Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    boost::ptr_map<int, GroupBase>::iterator itGroup1 = mGroupMap.find(rIdentGroup1);
    if (itGroup1 == mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupUnite] Group1 with the given identifier does not exist.");
    boost::ptr_map<int, GroupBase>::iterator itGroup2 = mGroupMap.find(rIdentGroup2);
    if (itGroup2 == mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupUnite] Group2 with the given identifier does not exist.");

    int groupNumber = GetUnusedId(mGroupMap);
    // insert the member
    mGroupMap.insert(groupNumber, (*itGroup1).second->Unite((*itGroup2).second));

    return groupNumber;
}

int NuTo::StructureBase::GroupDifference(int rIdentGroup1, int rIdentGroup2)
{
    NuTo::Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    boost::ptr_map<int, GroupBase>::iterator itGroup1 = mGroupMap.find(rIdentGroup1);
    if (itGroup1 == mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupDifference] Group1 with the given identifier does not exist.");
    boost::ptr_map<int, GroupBase>::iterator itGroup2 = mGroupMap.find(rIdentGroup2);
    if (itGroup2 == mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupDifference] Group2 with the given identifier does not exist.");

    int groupNumber = GetUnusedId(mGroupMap);

    // insert the member
    mGroupMap.insert(groupNumber, (*itGroup1).second->Difference((*itGroup2).second));
    return groupNumber;
}

int NuTo::StructureBase::GroupIntersection(int rIdentGroup1, int rIdentGroup2)
{
    NuTo::Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    boost::ptr_map<int, GroupBase>::iterator itGroup1 = mGroupMap.find(rIdentGroup1);
    if (itGroup1 == mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupIntersection] Group1 with the given identifier does not exist.");
    boost::ptr_map<int, GroupBase>::iterator itGroup2 = mGroupMap.find(rIdentGroup2);
    if (itGroup2 == mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupIntersection] Group2 with the given identifier does not exist.");

    int groupNumber = GetUnusedId(mGroupMap);
    // insert the member
    mGroupMap.insert(groupNumber, (*itGroup1).second->Intersection((*itGroup2).second));
    return groupNumber;
}

int NuTo::StructureBase::GroupSymmetricDifference(int rIdentGroup1, int rIdentGroup2)
{
    boost::ptr_map<int, GroupBase>::iterator itGroup1 = mGroupMap.find(rIdentGroup1);
    if (itGroup1 == mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupSymmetricDifference] Group1 with the given identifier does not exist.");
    boost::ptr_map<int, GroupBase>::iterator itGroup2 = mGroupMap.find(rIdentGroup2);
    if (itGroup2 == mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupSymmetricDifference] Group2 with the given identifier does not exist.");

    int groupNumber = GetUnusedId(mGroupMap);

    mGroupMap.insert(groupNumber, (*itGroup1).second->SymmetricDifference((*itGroup2).second));
    return groupNumber;
}

int NuTo::StructureBase::GroupGetNumMembers(int rIdentGroup) const
{
    NuTo::Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    boost::ptr_map<int, GroupBase>::const_iterator itGroup = mGroupMap.find(rIdentGroup);
    if (itGroup == mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupGetNumMembers] Group with the given identifier does not exist.");

    return (*itGroup).second->GetNumMembers();
}

std::vector<int> NuTo::StructureBase::GroupGetMemberIds(int rIdentGroup) const
{
    NuTo::Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    boost::ptr_map<int, GroupBase>::const_iterator itGroup = mGroupMap.find(rIdentGroup);
    if (itGroup == mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupGetNumMembers] Group with the given identifier does not exist.");

    return (*itGroup).second->GetMemberIds();
}

bool NuTo::StructureBase::GroupContainsMember(int rIdentGroup, int rMember) const
{
    NuTo::Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    boost::ptr_map<int, GroupBase>::const_iterator itGroup = mGroupMap.find(rIdentGroup);
    if (itGroup == mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupContains] Group with the given identifier does not exist.");

    return (*itGroup).second->Contain(rMember);
}


NuTo::GroupBase* NuTo::StructureBase::GroupGetGroupPtr(int rIdent)
{
    boost::ptr_map<int, GroupBase>::iterator it = this->mGroupMap.find(rIdent);
    if (it == this->mGroupMap.end())
    {
        throw NuTo::MechanicsException("[NuTo::StructureBase::GroupGetGroupPtr] Group does not exist.");
    }
    return it->second;
}


const NuTo::GroupBase* NuTo::StructureBase::GroupGetGroupPtr(int rIdent) const
{
    boost::ptr_map<int, GroupBase>::const_iterator it = this->mGroupMap.find(rIdent);
    if (it == this->mGroupMap.end())
    {
        throw NuTo::MechanicsException("[NuTo::StructureBase::GroupGetGroupPtr] Group does not exist.");
    }
    return it->second;
}
