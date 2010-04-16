// $Id$
#include "nuto/mechanics/groups/GroupBase.h"
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/loads/LoadNode.h"
#include "nuto/mechanics/loads/LoadNodeForces1D.h"
#include "nuto/mechanics/loads/LoadNodeForces3D.h"
#include "nuto/mechanics/loads/LoadNodeGroup.h"
#include "nuto/mechanics/loads/LoadNodeGroupForces1D.h"
#include "nuto/mechanics/loads/LoadNodeGroupForces3D.h"

// adds a force for a node
int NuTo::StructureBase::LoadCreateNodeForce(int rNodeIdent, const NuTo::FullMatrix<double>& rDirection, double rValue)
{
    // find node
    NodeBase* nodePtr;
    try
    {
        nodePtr = NodeGetNodePtr(rNodeIdent);
    }
    catch (NuTo::MechanicsException &e)
    {
        e.AddMessage("[NuTo::StructureBase::LoadCreateNodeForce] Node with the given identifier could not be found.");
        throw e;
    }
    catch (...)
    {
        throw MechanicsException("[NuTo::StructureBase::LoadCreateNodeForce] Node with the given identifier could not be found.");
    }

    return this->LoadCreateNodeForce(nodePtr,rDirection, rValue);
}

int NuTo::StructureBase::LoadCreateNodeForce(const NodeBase* rNode, const NuTo::FullMatrix<double>& rDirection, double rValue)
{
    //find unused integer id
    int id(0);
    boost::ptr_map<int,LoadBase>::iterator it = this->mLoadMap.find(id);
    while (it != this->mLoadMap.end())
    {
        id++;
        it = this->mLoadMap.find(id);
    }

    // create load
    LoadNode* loadPtr;
    switch (this->mDimension)
    {
    case 1:
        loadPtr = new NuTo::LoadNodeForces1D(rNode, rDirection(0,0), rValue);
        break;
    case 2:
        throw MechanicsException("[NuTo::StructureBase::LoadCreateNodeForce] To be implemented.");
        break;
    case 3:
        loadPtr = new NuTo::LoadNodeForces3D(rNode, rDirection, rValue);
        break;
    default:
        throw MechanicsException("[NuTo::StructureBase::LoadCreateNodeForce] Incorrect dimension of the structure.");
    }
    // insert load in load map
    this->mLoadMap.insert(id,loadPtr);
    return id;
}

// adds a force for a node group
int NuTo::StructureBase::LoadCreateNodeGroupForce(std::string rGroupIdent, const NuTo::FullMatrix<double>& rDirection, double rValue)
{
    // find group in map
    boost::ptr_map<std::string,GroupBase>::iterator itGroup = this->mGroupMap.find(rGroupIdent);
    if (itGroup == this->mGroupMap.end())
    {
        throw MechanicsException("[NuTo::Structure::LoadCreateNodeGroupForce] Group with the given identifier does not exist.");
    }

    // cast to node group (type check)
    Group<NodeBase> *nodeGroup = dynamic_cast<Group<NodeBase>*>(itGroup->second);
    if (nodeGroup == NULL)
    {
        throw MechanicsException("[NuTo::Structure::LoadCreateNodeGroupForce] Group is not a node group.");
    }

    return this->LoadCreateNodeGroupForce(nodeGroup,rDirection, rValue);
}

int NuTo::StructureBase::LoadCreateNodeGroupForce(const Group<NodeBase>* rNodeGroup, const NuTo::FullMatrix<double>& rDirection, double rValue)
{
    //find unused integer id
    int id(0);
    boost::ptr_map<int,LoadBase>::iterator it = this->mLoadMap.find(id);
    while (it != this->mLoadMap.end())
    {
        id++;
        it = this->mLoadMap.find(id);
    }

    // create load
    LoadNodeGroup* loadPtr;
    switch (this->mDimension)
    {
    case 1:
        loadPtr = new NuTo::LoadNodeGroupForces1D(rNodeGroup, rDirection(0,0), rValue);
        break;
    case 2:
        throw MechanicsException("[NuTo::StructureBase::LoadCreateNodeForce] To be implemented.");
        break;
    case 3:
        loadPtr = new NuTo::LoadNodeGroupForces3D(rNodeGroup, rDirection, rValue);
        break;
    default:
        throw MechanicsException("[NuTo::StructureBase::LoadCreateNodeForce] Incorrect dimension of the structure.");
    }
    // insert load in load map
    this->mLoadMap.insert(id,loadPtr);
    return id;
}

void NuTo::StructureBase::LoadDelete(int rIdent)
{
    // find load in map
    boost::ptr_map<int,LoadBase>::iterator it = this->mLoadMap.find(rIdent);
    if (it == this->mLoadMap.end())
    {
        throw NuTo::MechanicsException("[NuTo::StructureBase::Constitutive] Load does not exist.");
    }
    else
    {
        // delete load from map
        this->mLoadMap.erase(it);
    }
}

// get the pointer to a load from the load identifier
NuTo::LoadBase* NuTo::StructureBase::LoadGetLoadPtr(int rIdent)
{
    boost::ptr_map<int,LoadBase>::iterator it = this->mLoadMap.find(rIdent);
    if (it == this->mLoadMap.end())
    {
        throw NuTo::MechanicsException("[NuTo::StructureBase::ConstitutiveLaw] load does not exist.");
    }
    return it->second;

}

// get the pointer to a load from the load identifier
const NuTo::LoadBase* NuTo::StructureBase::LoadGetLoadPtr(int rIdent) const
{
    boost::ptr_map<int,LoadBase>::const_iterator it = this->mLoadMap.find(rIdent);
    if (it == this->mLoadMap.end())
    {
        throw NuTo::MechanicsException("[NuTo::StructureBase::ConstitutiveLaw] load does not exist.");
    }
    return it->second;

}
