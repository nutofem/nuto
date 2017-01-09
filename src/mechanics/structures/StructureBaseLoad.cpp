// $Id$
#include "mechanics/structures/StructureBase.h"
#include "mechanics/groups/Group.h"
#include "mechanics/integrationtypes/IntegrationTypeBase.h"
#include "mechanics/loads/LoadNode.h"
#include "mechanics/loads/LoadNodeForces1D.h"
#include "mechanics/loads/LoadNodeForces2D.h"
#include "mechanics/loads/LoadNodeForces3D.h"
#include "mechanics/loads/LoadNodeGroup.h"
#include "mechanics/loads/LoadNodeGroupForces1D.h"
#include "mechanics/loads/LoadNodeGroupForces2D.h"
#include "mechanics/loads/LoadNodeGroupForces3D.h"
#include "mechanics/loads/LoadNodeHeatFlux1D.h"
#include "mechanics/loads/LoadSurfaceConstDirection2D.h"
#include "mechanics/loads/LoadSurfaceConstDirection3D.h"
#include "mechanics/loads/LoadSurfacePressure2D.h"
#include "mechanics/loads/LoadSurfacePressureFunction2D.h"
#include "mechanics/loads/LoadSurfacePressure3D.h"

// adds a force for a node
int NuTo::StructureBase::LoadCreateNodeForce(int rLoadCase, int rNodeIdent, const Eigen::MatrixXd& rDirection, double rValue)
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
        throw;
    }
    catch (...)
    {
        throw MechanicsException("[NuTo::StructureBase::LoadCreateNodeForce] Node with the given identifier could not be found.");
    }

    return this->LoadCreateNodeForce(rLoadCase,nodePtr,rDirection, rValue);
}

int NuTo::StructureBase::LoadCreateNodeForce(int rLoadCase, const NodeBase* rNode, const Eigen::MatrixXd& rDirection, double rValue)
{
    if (rLoadCase>=mNumLoadCases)
    	throw MechanicsException("[NuTo::StructureBase::LoadCreateNodeForce] Load case number larger than total number of load cases. Use myStructure.SetNumLoadCases(num) to set the maximum number");
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

int NuTo::StructureBase::LoadCreateNodeHeatFlux(int rLoadCase, int rNodeIdent,
        const Eigen::MatrixXd& rDirection, double rValue)
{
    // find node
    NodeBase* nodePtr;
    try
    {
        nodePtr = NodeGetNodePtr(rNodeIdent);
    }
    catch (NuTo::MechanicsException &e)
    {
        e.AddMessage(__PRETTY_FUNCTION__, "Node with the given identifier could not be found.");
        throw;
    }
    catch (...)
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Node with the given identifier could not be found.");
    }

    return this->LoadCreateNodeHeatFlux(rLoadCase,nodePtr,rDirection, rValue);
}

int NuTo::StructureBase::LoadCreateNodeHeatFlux(int rLoadCase, const NodeBase* rNode,
        const Eigen::MatrixXd& rDirection, double rValue)
{
    if (rLoadCase>=mNumLoadCases)
        throw MechanicsException(__PRETTY_FUNCTION__, "Load case number larger than total number of load cases. Use myStructure.SetNumLoadCases(num) to set the maximum number");
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
        loadPtr = new NuTo::LoadNodeHeatFlux1D(rLoadCase, rNode, rDirection(0,0), rValue);
        break;
    case 2:
        //loadPtr = new NuTo::LoadNodeForces2D(rLoadCase, rNode, rDirection, rValue);
        throw MechanicsException(__PRETTY_FUNCTION__, "Boundary heat flux for 2D not yet implemented.");
        break;
    case 3:
        //loadPtr = new NuTo::LoadNodeForces3D(rLoadCase, rNode, rDirection, rValue);
        throw MechanicsException(__PRETTY_FUNCTION__, "Boundary heat flux for 3D not yet implemented.");
        break;
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Incorrect dimension of the structure.");
    }
    // insert load in load map
    this->mLoadMap.insert(id,loadPtr);
    return id;
}

// adds a force for a node group
int NuTo::StructureBase::LoadCreateNodeGroupForce(int rLoadCase, int rGroupIdent, const Eigen::MatrixXd& rDirection, double rValue)
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

int NuTo::StructureBase::LoadCreateNodeGroupForce(int rLoadCase, const Group<NodeBase>* rNodeGroup, const Eigen::MatrixXd& rDirection, double rValue)
{
    if (rLoadCase>=mNumLoadCases)
    	throw MechanicsException("[NuTo::StructureBase::LoadCreateNodeGroupForce] Load case number larger than total number of load cases. Use myStructure.SetNumLoadCases(num) to set the maximum number");
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

int NuTo::StructureBase::LoadSurfaceConstDirectionCreate3D(int rLoadCase, int rElementGroupId, int rNodeGroupId, const Eigen::VectorXd& rLoadVector)
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

int NuTo::StructureBase::LoadSurfaceConstDirectionCreate2D(int rLoadCase, int rElementGroupId, int rNodeGroupId, const Eigen::VectorXd& rLoadVector)
{
    if (rLoadCase>=mNumLoadCases)
    	throw MechanicsException("[NuTo::StructureBase::LoadSurfaceConstDirectionCreate2D] Load case number larger than total number of load cases. Use myStructure.SetNumLoadCases(num) to set the maximum number");
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
    if (rLoadCase>=mNumLoadCases)
    	throw MechanicsException("[NuTo::StructureBase::LoadSurfacePressureCreate3D] Load case number larger than total number of load cases. Use myStructure.SetNumLoadCases(num) to set the maximum number");
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
    if (rLoadCase>=mNumLoadCases)
    	throw MechanicsException("[NuTo::StructureBase::LoadSurfacePressureCreate2D] Load case number larger than total number of load cases. Use myStructure.SetNumLoadCases(num) to set the maximum number");
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

int NuTo::StructureBase::LoadSurfacePressureFunctionCreate2D(int rLoadCase,
                                                             int rElementGroupId,
                                                             int rNodeGroupId,
                                                             const std::function<NuTo::FullVector<double,2>(NuTo::FullVector<double,2>)> &rLoadFunction)
{
    if (rLoadCase>=mNumLoadCases)
        throw MechanicsException("[NuTo::StructureBase::LoadSurfacePressureFunctionCreate2D] Load case number larger than total number of load cases. Use myStructure.SetNumLoadCases(num) to set the maximum number");
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
    loadPtr = new NuTo::LoadSurfacePressureFunction2D(rLoadCase, &(*this), rElementGroupId, rNodeGroupId, rLoadFunction);

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
