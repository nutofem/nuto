#include "mechanics/structures/StructureBase.h"
#include "mechanics/groups/Group.h"
#include "mechanics/integrationtypes/IntegrationTypeBase.h"
#include "mechanics/loads/LoadNode.h"
#include "mechanics/loads/LoadNodeScalarSource.h"
#include "mechanics/loads/LoadNodeForces1D.h"
#include "mechanics/loads/LoadNodeForces2D.h"
#include "mechanics/loads/LoadNodeForces3D.h"
#include "mechanics/loads/LoadNodeGroup.h"
#include "mechanics/loads/LoadNodeGroupForces1D.h"
#include "mechanics/loads/LoadNodeGroupForces2D.h"
#include "mechanics/loads/LoadNodeGroupForces3D.h"
#include "mechanics/loads/LoadSurfaceConstDirection2D.h"
#include "mechanics/loads/LoadSurfaceConstDirection3D.h"
#include "mechanics/loads/LoadSurfacePressure2D.h"
#include "mechanics/loads/LoadSurfacePressureFunction2D.h"
#include "mechanics/loads/LoadSurfacePressure3D.h"

using namespace NuTo;

int StructureBase::LoadCreateScalarSource(int rNodeIdent, double rValue)
{
    // find node
    NodeBase* nodePtr;
    nodePtr = NodeGetNodePtr(rNodeIdent);
    return this->LoadCreateScalarSource(nodePtr, rValue);
}


int StructureBase::LoadCreateScalarSource(const NodeBase* rNode, double rValue)
{
    int id = GetUnusedId(mLoadMap);

    // create load
    LoadNode* loadPtr;
    loadPtr = new LoadNodeScalarSource(rNode, rValue);

    // insert load in load map
    this->mLoadMap.insert(id,loadPtr);
    return id;
}

int StructureBase::LoadCreateNodeForce(int rNodeIdent, const Eigen::MatrixXd& rDirection, double rValue)
{
    auto nodePtr = NodeGetNodePtr(rNodeIdent);
    return LoadCreateNodeForce(nodePtr, rDirection, rValue);
}


int StructureBase::LoadCreateNodeForce(const NodeBase* rNode, const Eigen::MatrixXd& rDirection, double rValue)
{
    int id = GetUnusedId(mLoadMap);

    // create load
    LoadNode* loadPtr;
    switch (mDimension)
    {
    case 1:
        loadPtr = new LoadNodeForces1D(rNode, rDirection(0, 0), rValue);
        break;
    case 2:
        loadPtr = new LoadNodeForces2D(rNode, rDirection, rValue);
        break;
    case 3:
        loadPtr = new LoadNodeForces3D(rNode, rDirection, rValue);
        break;
    default:
        throw Exception(__PRETTY_FUNCTION__, "Incorrect dimension of the structure.");
    }
    // insert load in load map
    mLoadMap.insert(id, loadPtr);
    return id;
}


int StructureBase::LoadCreateNodeGroupForce(int rGroupIdent, const Eigen::MatrixXd& rDirection, double rValue)
{
    // find group in map
    auto itGroup = mGroupMap.find(rGroupIdent);
    if (itGroup == mGroupMap.end())
    {
        throw Exception(__PRETTY_FUNCTION__, "Group with the given identifier does not exist.");
    }

    // cast to node group (type check)
    Group<NodeBase>* nodeGroup = dynamic_cast<Group<NodeBase>*>(itGroup->second);
    if (nodeGroup == nullptr)
    {
        throw Exception(__PRETTY_FUNCTION__, "Group is not a node group.");
    }

    return LoadCreateNodeGroupForce(nodeGroup, rDirection, rValue);
}


int StructureBase::LoadCreateNodeGroupForce(const Group<NodeBase>* rNodeGroup, const Eigen::MatrixXd& rDirection,
                                            double rValue)
{
    int id = GetUnusedId(mLoadMap);

    // create load
    LoadNodeGroup* loadPtr;
    switch (mDimension)
    {
    case 1:
        loadPtr = new LoadNodeGroupForces1D(rNodeGroup, rDirection(0, 0), rValue);
        break;
    case 2:
        loadPtr = new LoadNodeGroupForces2D(rNodeGroup, rDirection, rValue);
        break;
    case 3:
        loadPtr = new LoadNodeGroupForces3D(rNodeGroup, rDirection, rValue);
        break;
    default:
        throw Exception(__PRETTY_FUNCTION__, "Incorrect dimension of the structure.");
    }
    // insert load in load map
    mLoadMap.insert(id, loadPtr);

    return id;
}


int StructureBase::LoadSurfaceConstDirectionCreate3D(int rElementGroupId, int rNodeGroupId,
                                                     const Eigen::VectorXd& rLoadVector)
{
    int id = GetUnusedId(mLoadMap);

    // create load
    LoadSurfaceBase3D* loadPtr;
    loadPtr = new LoadSurfaceConstDirection3D(this, rElementGroupId, rNodeGroupId, rLoadVector);

    // insert load in load map
    mLoadMap.insert(id, loadPtr);
    return id;
}


int StructureBase::LoadSurfaceConstDirectionCreate2D(int rElementGroupId, int rNodeGroupId,
                                                     const Eigen::VectorXd& rLoadVector)
{
    int id = GetUnusedId(mLoadMap);

    // create load
    LoadSurfaceBase2D* loadPtr;
    loadPtr = new LoadSurfaceConstDirection2D(this, rElementGroupId, rNodeGroupId, rLoadVector);

    // insert load in load map
    mLoadMap.insert(id, loadPtr);
    return id;
}


int StructureBase::LoadSurfacePressureCreate3D(int rElementGroupId, int rNodeGroupId, double rPressure)
{
    int id = GetUnusedId(mLoadMap);

    // create load
    LoadSurfaceBase3D* loadPtr;
    loadPtr = new LoadSurfacePressure3D(this, rElementGroupId, rNodeGroupId, rPressure);

    // insert load in load map
    mLoadMap.insert(id, loadPtr);
    return id;
}


int StructureBase::LoadSurfacePressureCreate2D(int rElementGroupId, int rNodeGroupId, double rPressure)
{
    int id = GetUnusedId(mLoadMap);

    // create load
    LoadSurfaceBase2D* loadPtr;
    loadPtr = new LoadSurfacePressure2D(this, rElementGroupId, rNodeGroupId, rPressure);

    // insert load in load map
    mLoadMap.insert(id, loadPtr);
    return id;
}


int StructureBase::LoadSurfacePressureFunctionCreate2D(
        int rElementGroupId, int rNodeGroupId, const std::function<Eigen::Vector2d(Eigen::Vector2d)>& rLoadFunction)
{
    int id = GetUnusedId(mLoadMap);

    // create load
    LoadSurfaceBase2D* loadPtr;
    loadPtr = new LoadSurfacePressureFunction2D(this, rElementGroupId, rNodeGroupId, rLoadFunction);

    // insert load in load map
    mLoadMap.insert(id, loadPtr);
    return id;
}
