// $Id$
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"
#include "nuto/mechanics/loads/LoadNode.h"
#include "nuto/mechanics/loads/LoadNodeForces1D.h"
#include "nuto/mechanics/loads/LoadNodeForces2D.h"
#include "nuto/mechanics/loads/LoadNodeForces3D.h"
#include "nuto/mechanics/loads/LoadNodeGroup.h"
#include "nuto/mechanics/loads/LoadNodeGroupForces1D.h"
#include "nuto/mechanics/loads/LoadNodeGroupForces2D.h"
#include "nuto/mechanics/loads/LoadNodeGroupForces3D.h"
#include "nuto/mechanics/loads/LoadSurfaceConstDirection2D.h"
#include "nuto/mechanics/loads/LoadSurfaceConstDirection3D.h"
#include "nuto/mechanics/loads/LoadSurfacePressure2D.h"
#include "nuto/mechanics/loads/LoadSurfacePressure3D.h"

// adds a force for a node
int NuTo::StructureBase::LoadCreateNodeForce(int rLoadCase, int rNodeIdent, const NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rDirection, double rValue)
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

    return this->LoadCreateNodeForce(rLoadCase,nodePtr,rDirection, rValue);
}

int NuTo::StructureBase::LoadCreateNodeForce(int rLoadCase, const NodeBase* rNode, const NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rDirection, double rValue)
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
        loadPtr = new NuTo::LoadNodeForces1D(rLoadCase, rNode, rDirection(0,0), rValue);
        break;
    case 2:
        loadPtr = new NuTo::LoadNodeForces2D(rLoadCase, rNode, rDirection, rValue);
        break;
    case 3:
        loadPtr = new NuTo::LoadNodeForces3D(rLoadCase, rNode, rDirection, rValue);
        break;
    default:
        throw MechanicsException("[NuTo::StructureBase::LoadCreateNodeForce] Incorrect dimension of the structure.");
    }
    // insert load in load map
    this->mLoadMap.insert(id,loadPtr);
    return id;
}

// adds a force for a node group
int NuTo::StructureBase::LoadCreateNodeGroupForce(int rLoadCase, int rGroupIdent, const NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rDirection, double rValue)
{
    // find group in map
    boost::ptr_map<int,GroupBase>::iterator itGroup = this->mGroupMap.find(rGroupIdent);
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

    return this->LoadCreateNodeGroupForce(rLoadCase, nodeGroup,rDirection, rValue);
}

int NuTo::StructureBase::LoadCreateNodeGroupForce(int rLoadCase, const Group<NodeBase>* rNodeGroup, const NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rDirection, double rValue)
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
        loadPtr = new NuTo::LoadNodeGroupForces1D(rLoadCase, rNodeGroup, rDirection(0,0), rValue);
        break;
    case 2:
        loadPtr = new NuTo::LoadNodeGroupForces2D(rLoadCase, rNodeGroup, rDirection, rValue);
        break;
    case 3:
        loadPtr = new NuTo::LoadNodeGroupForces3D(rLoadCase, rNodeGroup, rDirection, rValue);
        break;
    default:
        throw MechanicsException("[NuTo::StructureBase::LoadCreateNodeForce] Incorrect dimension of the structure.");
    }
    // insert load in load map
    this->mLoadMap.insert(id,loadPtr);
    return id;
}

int NuTo::StructureBase::LoadSurfaceConstDirectionCreate3D(int rLoadCase, int rElementGroupId, int rNodeGroupId, const NuTo::FullVector<double,Eigen::Dynamic>& rLoadVector)
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
    LoadSurfaceBase3D* loadPtr;
    loadPtr = new NuTo::LoadSurfaceConstDirection3D(rLoadCase, &(*this), rElementGroupId, rNodeGroupId, rLoadVector);

    // insert load in load map
    this->mLoadMap.insert(id,loadPtr);
    return id;
}

int NuTo::StructureBase::LoadSurfaceConstDirectionCreate2D(int rLoadCase, int rElementGroupId, int rNodeGroupId, const NuTo::FullVector<double,Eigen::Dynamic>& rLoadVector)
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
    LoadSurfaceBase2D* loadPtr;
    loadPtr = new NuTo::LoadSurfaceConstDirection2D(rLoadCase, &(*this), rElementGroupId, rNodeGroupId, rLoadVector);

    // insert load in load map
    this->mLoadMap.insert(id,loadPtr);
    return id;
}

int NuTo::StructureBase::LoadSurfacePressureCreate3D(int rLoadCase, int rElementGroupId, int rNodeGroupId, double rPressure)
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
    LoadSurfaceBase3D* loadPtr;
    loadPtr = new NuTo::LoadSurfacePressure3D(rLoadCase, &(*this), rElementGroupId, rNodeGroupId, rPressure);

    // insert load in load map
    this->mLoadMap.insert(id,loadPtr);
    return id;
}

int NuTo::StructureBase::LoadSurfacePressureCreate2D(int rLoadCase, int rElementGroupId, int rNodeGroupId, double rPressure)
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
    LoadSurfaceBase2D* loadPtr;
    loadPtr = new NuTo::LoadSurfacePressure2D(rLoadCase, &(*this), rElementGroupId, rNodeGroupId, rPressure);

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
