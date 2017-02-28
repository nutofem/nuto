// $Id$

#include <boost/assign/ptr_map_inserter.hpp>

#include "base/Timer.h"


#include "mechanics/elements/ElementOutputFullMatrixDouble.h"
#include "mechanics/elements/ElementOutputFullVectorDouble.h"
#include "mechanics/elements/ElementOutputBlockVectorDouble.h"
#include "mechanics/elements/ElementOutputBlockVectorInt.h"
#include "mechanics/elements/ElementOutputVectorInt.h"

#include "mechanics/structures/StructureBase.h"
#include "mechanics/structures/StructureOutputBlockVector.h"
#include "mechanics/elements/ElementBase.h"
#include "mechanics/elements/ElementEnum.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/groups/Group.h"
#include "mechanics/groups/GroupEnum.h"

void NuTo::StructureBase::NodeSetDisplacements(int rNode, const Eigen::VectorXd& rDisplacements)
{
    NuTo::Timer(__FUNCTION__, GetShowTime(), GetLogger());
    NodeBase* nodePtr=NodeGetNodePtr(rNode);
    this->mUpdateTmpStaticDataRequired=true;

    if (rDisplacements.cols()!=1)
        throw MechanicsException(__PRETTY_FUNCTION__, "Displacement matrix has to have a single column.");



    try
    {
        if (rDisplacements.rows() <= 0 or rDisplacements.rows() > 3)
            throw MechanicsException(__PRETTY_FUNCTION__, "The number of displacement components is either 1, 2 or 3.");

        nodePtr->Set(Node::eDof::DISPLACEMENTS, rDisplacements);

    }
    catch(NuTo::MechanicsException & b)
    {
        b.AddMessage("[NuTo::StructureBase::NodeSetDisplacements] Error setting displacements.");
        throw;
    }
    catch(...)
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Error setting displacements of node (unspecified exception).");
    }
}


NuTo::StructureOutputBlockVector NuTo::StructureBase::NodeExtractDofValues() const
{
    return NodeExtractDofValues(0);
}

void NuTo::StructureBase::NodeSetDisplacements(int rNode, int rTimeDerivative, const Eigen::VectorXd& rDisplacements)
{
    NuTo::Timer(__FUNCTION__, GetShowTime(), GetLogger());
    NodeBase* nodePtr=NodeGetNodePtr(rNode);
    this->mUpdateTmpStaticDataRequired=true;

    if (rDisplacements.cols()!=1)
        throw MechanicsException(__PRETTY_FUNCTION__, "Displacement matrix has to have a single column.");
    if (nodePtr->GetNumTimeDerivatives(Node::eDof::DISPLACEMENTS)<rTimeDerivative)
        throw MechanicsException(__PRETTY_FUNCTION__, "number of time derivatives stored at node is less than the required value.");
    try
    {
        if (rDisplacements.rows() <= 0 or rDisplacements.rows() > 3)
            throw MechanicsException(__PRETTY_FUNCTION__, "The number of displacement components is either 1, 2 or 3.");

        nodePtr->Set(Node::eDof::DISPLACEMENTS, rTimeDerivative, rDisplacements);
    }
    catch(NuTo::MechanicsException & b)
    {
        b.AddMessage("[NuTo::StructureBase::NodeSetDisplacements] Error setting displacements.");
        throw;
    }
    catch(...)
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Error setting displacements of node (unspecified exception).");
    }
}

void NuTo::StructureBase::NodeSetRotations(int rNode, const Eigen::VectorXd& rRotations)
{
    NuTo::Timer(__FUNCTION__, GetShowTime(), GetLogger());
    NodeBase* nodePtr=NodeGetNodePtr(rNode);
    this->mUpdateTmpStaticDataRequired=true;

    if (rRotations.cols()!=1)
        throw MechanicsException(__PRETTY_FUNCTION__, "rotation matrix has to have a single column.");
    try
    {
        if (rRotations.rows() != 1 and rRotations.rows() != 3)
            throw MechanicsException(__PRETTY_FUNCTION__, "The number of rotation components is either 1, 3.");

        nodePtr->Set(Node::eDof::ROTATIONS, rRotations);
    }
    catch(NuTo::MechanicsException & b)
    {
        b.AddMessage("[NuTo::StructureBase::NodeSetRotations] Error setting rotations.");
        throw;
    }
    catch(...)
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Error setting rotations of node (unspecified exception).");
    }
}

void NuTo::StructureBase::NodeGroupSetDisplacements(int rGroupIdent, const Eigen::VectorXd& rDisplacements)
{
    NuTo::Timer(__FUNCTION__, GetShowTime(), GetLogger());
    this->mUpdateTmpStaticDataRequired=true;
    if (rDisplacements.cols()!=1)
        throw MechanicsException(__PRETTY_FUNCTION__, "Displacement matrix has to have a single column.");

    boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rGroupIdent);
    if (itGroup==mGroupMap.end())
        throw MechanicsException(__PRETTY_FUNCTION__, "Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::eGroupId::Nodes)
        throw MechanicsException(__PRETTY_FUNCTION__, "Group is not a node group.");
    Group<NodeBase> *nodeGroup = dynamic_cast<Group<NodeBase>*>(itGroup->second);
    assert(nodeGroup!=0);

    for (Group<NodeBase>::iterator itNode=nodeGroup->begin(); itNode!=nodeGroup->end();itNode++)
    {
        try
        {
            if (rDisplacements.rows() <= 0 or rDisplacements.rows() > 3)
                throw MechanicsException(__PRETTY_FUNCTION__, "The number of displacement components is either 1, 2 or 3.");

            itNode->second->Set(Node::eDof::DISPLACEMENTS, rDisplacements);
        }
        catch(NuTo::MechanicsException & b)
        {
            b.AddMessage("[NuTo::StructureBase::NodeGroupSetDisplacements] Error setting displacements.");
            throw;
        }
        catch(...)
        {
            throw MechanicsException(__PRETTY_FUNCTION__, "Error setting displacements of node (unspecified exception).");
        }
    }
}

void NuTo::StructureBase::NodeGroupSetDisplacements(int rGroupIdent, int rTimeDerivative, const Eigen::VectorXd& rDisplacements)
{
    NuTo::Timer(__FUNCTION__, GetShowTime(), GetLogger());
    this->mUpdateTmpStaticDataRequired=true;
    if (rDisplacements.cols()!=1)
        throw MechanicsException(__PRETTY_FUNCTION__, "Displacement matrix has to have a single column.");

    boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rGroupIdent);
    if (itGroup==mGroupMap.end())
        throw MechanicsException(__PRETTY_FUNCTION__, "Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::eGroupId::Nodes)
        throw MechanicsException(__PRETTY_FUNCTION__, "Group is not a node group.");
    Group<NodeBase> *nodeGroup = dynamic_cast<Group<NodeBase>*>(itGroup->second);
    assert(nodeGroup!=0);

    for (Group<NodeBase>::iterator itNode=nodeGroup->begin(); itNode!=nodeGroup->end();itNode++)
    {
        try
        {
            if (itNode->second->GetNumTimeDerivatives(Node::eDof::DISPLACEMENTS)<rTimeDerivative)
                throw MechanicsException(__PRETTY_FUNCTION__, "does not have a sufficient number of time derivatives.");

            if (rDisplacements.rows() <= 0 or rDisplacements.rows() > 3)
                throw MechanicsException(__PRETTY_FUNCTION__, "The number of displacement components is either 1, 2 or 3.");

            itNode->second->Set(Node::eDof::DISPLACEMENTS, rTimeDerivative, rDisplacements);
        }
        catch(NuTo::MechanicsException & b)
        {
            b.AddMessage("[NuTo::StructureBase::NodeGroupSetDisplacements] Error setting displacements.");
            throw;
        }
        catch(...)
        {
            throw MechanicsException(__PRETTY_FUNCTION__, "Error setting displacements of node (unspecified exception).");
        }
    }
}

void NuTo::StructureBase::NodeSetTemperature(int rNode, double rTemperature)
{
    NodeBase* nodePtr = NodeGetNodePtr(rNode);
    this->mUpdateTmpStaticDataRequired = true;
    nodePtr->Set(Node::eDof::TEMPERATURE, rTemperature);
}

void NuTo::StructureBase::NodeSetTemperature(int rNode, int rTimeDerivative, double rTemperature)
{
    NodeBase* nodePtr = NodeGetNodePtr(rNode);
    this->mUpdateTmpStaticDataRequired = true;
    if (nodePtr->GetNumTimeDerivatives(Node::eDof::TEMPERATURE)<rTimeDerivative)
        throw MechanicsException(__PRETTY_FUNCTION__,
                "Number of time derivatives stored at node is less than the required value.");
    nodePtr->Set(Node::eDof::TEMPERATURE, rTimeDerivative, rTemperature);
}

void NuTo::StructureBase::NodeGroupGetMembers(int rGroupId, std::vector<int>& rMembers)
{
    NuTo::Timer(__FUNCTION__, GetShowTime(), GetLogger());

    boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rGroupId);
    if (itGroup==mGroupMap.end())
        throw MechanicsException(__PRETTY_FUNCTION__, "Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::eGroupId::Nodes)
        throw MechanicsException(__PRETTY_FUNCTION__, "Group is not an node group.");
    Group<NodeBase> *nodeGroup = itGroup->second->AsGroupNode();
    assert(nodeGroup!=0);

    rMembers.resize(nodeGroup->GetNumMembers());
    int countNode(0);
    for (Group<NodeBase>::const_iterator itNode=nodeGroup->begin(); itNode!=nodeGroup->end();itNode++,countNode++)
    {
        rMembers[countNode] = itNode->first;
    }
}


void NuTo::StructureBase::NodeGetDisplacements(int rNode, Eigen::VectorXd& rDisplacements)const
{
    this->NodeGetDisplacements(rNode,0,rDisplacements);
}


void NuTo::StructureBase::NodeGetDisplacements(int rNode, int rTimeDerivative, Eigen::VectorXd& rDisplacements)const
{
    NuTo::Timer(__FUNCTION__, GetShowTime(), GetLogger());

    const NodeBase* nodePtr = NodeGetNodePtr(rNode);

    try
    {
        if (nodePtr->GetNum(Node::eDof::DISPLACEMENTS) == 0)
            throw MechanicsException(__PRETTY_FUNCTION__, "Node has no displacements.");

        rDisplacements = nodePtr->Get(Node::eDof::DISPLACEMENTS, rTimeDerivative);
    }
    catch(NuTo::MechanicsException & b)
    {
        b.AddMessage("[NuTo::StructureBase::NodeGetDisplacements] Error getting displacements.");
        throw;
    }
    catch(...)
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Error getting displacements of node (unspecified exception).");
    }
}


std::vector<int> NuTo::StructureBase::NodeGetDofIds(const int rNodeId, NuTo::Node::eDof rDof)const
{
    NuTo::Timer(__FUNCTION__, GetShowTime(), GetLogger());

    const NodeBase* nodePtr = NodeGetNodePtr(rNodeId);

    try
    {
        int numDofIds = nodePtr->GetNum(rDof);

        if (numDofIds == 0)
            throw MechanicsException(__PRETTY_FUNCTION__, "Node does not have the requested dof.");

        std::vector<int> dofIds(numDofIds);
        for (int i = 0; i < numDofIds; ++i)
            dofIds[i] = nodePtr->GetDof(rDof, i);

        return dofIds;

    }
    catch(NuTo::MechanicsException & b)
    {
        b.AddMessage(__PRETTY_FUNCTION__, "Error getting the requested dof identifiers.");
        throw;
    }
    catch(...)
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Error getting the requested dof identifiers.");
    }
}

void NuTo::StructureBase::NodeGetRotations(int rNode, Eigen::VectorXd& rRotations)const
{
    NuTo::Timer(__FUNCTION__, GetShowTime(), GetLogger());

    const NodeBase* nodePtr = NodeGetNodePtr(rNode);

    try
    {
        if (nodePtr->GetNum(Node::eDof::ROTATIONS) != 1 and nodePtr->GetNum(Node::eDof::ROTATIONS) != 3)
            throw MechanicsException(__PRETTY_FUNCTION__, "Node has neither 1(2D) or 3(3D) rotations.");

        rRotations = nodePtr->Get(Node::eDof::ROTATIONS);
    }
    catch(NuTo::MechanicsException & b)
    {
        b.AddMessage("[NuTo::StructureBase::NodeGetRotations] Error getting rotations.");
        throw;
    }
    catch(...)
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Error getting rotations of node (unspecified exception).");
    }
}
void NuTo::StructureBase::NodeGroupGetDisplacements(int rGroupIdent, Eigen::MatrixXd& rDisplacements)
{
    NuTo::Timer(__FUNCTION__, GetShowTime(), GetLogger());

    boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rGroupIdent);
    if (itGroup==mGroupMap.end())
        throw MechanicsException(__PRETTY_FUNCTION__, "Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::eGroupId::Nodes)
        throw MechanicsException(__PRETTY_FUNCTION__, "Group is not a node group.");
    Group<NodeBase> *nodeGroup = itGroup->second->AsGroupNode();
    assert(nodeGroup!=0);

    //all nodes have to have the same dimension
    if(nodeGroup->GetNumMembers()<1)
        throw MechanicsException(__PRETTY_FUNCTION__, "Group has no members.");

    int numDisp= nodeGroup->begin()->second->GetNum(Node::eDof::DISPLACEMENTS);
    //resize the matrix
    rDisplacements.resize(nodeGroup->GetNumMembers(),numDisp);

    int theNode(0);
    for (Group<NodeBase>::iterator itNode=nodeGroup->begin(); itNode!=nodeGroup->end();itNode++, theNode++)
    {
        try
        {
            if (numDisp != 1 and numDisp != 2 and numDisp != 3)
                throw MechanicsException(__PRETTY_FUNCTION__, "The number of displacement components is either 1, 2 or 3.");

            rDisplacements.row(theNode) = itNode->second->Get(Node::eDof::DISPLACEMENTS).transpose();

        }
        catch(NuTo::MechanicsException & b)
        {
            b.AddMessage("[NuTo::StructureBase::NodeGroupGetDisplacements] Error getting displacements.");
            throw;
        }
        catch(...)
        {
            throw MechanicsException(__PRETTY_FUNCTION__, "Error getting displacements of node (unspecified exception).");
        }
    }
}

double NuTo::StructureBase::NodeGetTemperature(int rNode) const
{
    return this->NodeGetTemperature(rNode, 0);
}

double NuTo::StructureBase::NodeGetTemperature(int rNode, int rTimeDerivative) const
{
    const NodeBase* nodePtr = NodeGetNodePtr(rNode);
    if (nodePtr->GetNum(Node::eDof::TEMPERATURE) == 0)
        throw MechanicsException(__PRETTY_FUNCTION__, "Node doesn't have a temperature.");
    return nodePtr->Get(Node::eDof::TEMPERATURE, rTimeDerivative)[0];
}

void NuTo::StructureBase::NodeGetCoordinates(int rNode, Eigen::VectorXd& rCoordinates)const
{
    NuTo::Timer(__FUNCTION__, GetShowTime(), GetLogger());

    const NodeBase* nodePtr = NodeGetNodePtr(rNode);

    try
    {
        if (nodePtr->GetNum(Node::eDof::COORDINATES) == 0)
            throw MechanicsException(__PRETTY_FUNCTION__, "Node has no coordinates.");

        rCoordinates = nodePtr->Get(Node::eDof::COORDINATES);

    }
    catch(NuTo::MechanicsException & b)
    {
        b.AddMessage("[NuTo::StructureBase::NodeGetCoordinates] Error getting coordinates.");
        throw;
    }
    catch(...)
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Error getting coordinates of node (unspecified exception).");
    }
}

void NuTo::StructureBase::NodeGroupGetCoordinates(int rGroupIdent, Eigen::MatrixXd& rCoordinates)
{
    NuTo::Timer(__FUNCTION__, GetShowTime(), GetLogger());

    boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rGroupIdent);
    if (itGroup==mGroupMap.end())
        throw MechanicsException(__PRETTY_FUNCTION__, "Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::eGroupId::Nodes)
        throw MechanicsException(__PRETTY_FUNCTION__, "Group is not a node group.");
    Group<NodeBase> *nodeGroup = itGroup->second->AsGroupNode();
    assert(nodeGroup!=0);

    //all nodes have to have the same dimension
    if(nodeGroup->GetNumMembers()<1)
        throw MechanicsException(__PRETTY_FUNCTION__, "Group has no members.");

    int numCoords= nodeGroup->begin()->second->GetNum(Node::eDof::COORDINATES);
    //resize the matrix
    rCoordinates.resize(nodeGroup->GetNumMembers(),numCoords);
    int theNode(0);
    for (Group<NodeBase>::iterator itNode=nodeGroup->begin(); itNode!=nodeGroup->end();itNode++, theNode++)
    {
        try
        {
            if (numCoords != 1 and numCoords != 2 and numCoords != 3)
                throw MechanicsException(__PRETTY_FUNCTION__, "The number of coordinates components is either 1, 2 or 3.");

            rCoordinates.row(theNode) = itNode->second->Get(Node::eDof::COORDINATES).transpose();

        }
        catch(NuTo::MechanicsException & b)
        {
            b.AddMessage("[NuTo::StructureBase::NodeGroupGetCoordinates] Error getting coordinates.");
            throw;
        }
        catch(...)
        {
            throw MechanicsException(__PRETTY_FUNCTION__, "Error getting coordinates of node (unspecified exception).");
        }
    }
}

void NuTo::StructureBase::NodeGetNonlocalEqPlasticStrain(int rNode, Eigen::VectorXd& rNonlocalEqPlasticStrain)const
{
    NuTo::Timer(__FUNCTION__, GetShowTime(), GetLogger());
    const NodeBase* nodePtr = NodeGetNodePtr(rNode);

    try
    {
        if (nodePtr->GetNum(Node::eDof::NONLOCALEQPLASTICSTRAIN) != 2)
        {
            throw MechanicsException(__PRETTY_FUNCTION__, "Node does not have nonlocal equivalent plastic strains.");
        }
        rNonlocalEqPlasticStrain = nodePtr->Get(Node::eDof::NONLOCALEQPLASTICSTRAIN);
    }
    catch(NuTo::MechanicsException & b)
    {
        b.AddMessage("[NuTo::StructureBase::NodeGetNonlocalEqPlasticStrain] Error getting global damage.");
        throw;
    }
    catch(...)
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Error getting NodeGetNonlocalEqPlasticStrain of node (unspecified exception).");
    }
}

void NuTo::StructureBase::NodeGetNonlocalTotalStrain(int rNode, Eigen::VectorXd& rNonlocalTotalStrain)const
{
    NuTo::Timer(__FUNCTION__, GetShowTime(), GetLogger());

    const NodeBase* nodePtr = NodeGetNodePtr(rNode);

    try
    {
        int num = nodePtr->GetNum(Node::eDof::NONLOCALTOTALSTRAIN);
        if (num != 1 and num != 3 and num != 6)
            throw MechanicsException(__PRETTY_FUNCTION__, "Number of nonlocal total strain components is either 1, 3 or 6 .");

        rNonlocalTotalStrain = nodePtr->Get(Node::eDof::NONLOCALTOTALSTRAIN);

    }
    catch(NuTo::MechanicsException & b)
    {
        b.AddMessage("[NuTo::StructureBase::NodeGetNonlocalTotalStrain] Error getting nonlocal total strain.");
        throw;
    }
    catch(...)
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Error getting nonlocal total strain of node (unspecified exception).");
    }
}


void NuTo::StructureBase::NodeInternalForce(int rId, Eigen::VectorXd& rNodeForce)
{
    NuTo::Timer(__FUNCTION__, GetShowTime(), GetLogger());

    try
    {
        const NodeBase* nodePtr = NodeGetNodePtr(rId);
        NodeInternalForce(nodePtr,rNodeForce);
    }
    catch(NuTo::MechanicsException & b)
    {
        b.AddMessage("[NuTo::StructureBase::NodeGradientInternalPotential] Error getting gradient of internal potential.");
        throw;
    }
    catch(...)
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Error getting gradient of internal potential (unspecified exception).");
    }
}

void NuTo::StructureBase::NodeGroupInternalForce(int rGroupIdent, Eigen::VectorXd& rNodeForce)
{
    NuTo::Timer(__FUNCTION__, GetShowTime(), GetLogger());

    boost::ptr_map<int,GroupBase>::const_iterator itGroup = mGroupMap.find(rGroupIdent);
    if (itGroup==mGroupMap.end())
        throw MechanicsException(__PRETTY_FUNCTION__, "Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::eGroupId::Nodes)
        throw MechanicsException(__PRETTY_FUNCTION__, "Group is not a node group.");
    const Group<NodeBase> *nodeGroup = dynamic_cast<const Group<NodeBase>*>(itGroup->second);
    assert(nodeGroup!=0);

    Eigen::VectorXd nodeForceLocal;

    if (nodeGroup->GetNumMembers()==0)
        throw MechanicsException(__PRETTY_FUNCTION__, "Node group is empty.");
    rNodeForce.resize(nodeGroup->begin()->second->GetNum(Node::eDof::DISPLACEMENTS));
    rNodeForce.setZero();

    for (Group<NodeBase>::const_iterator itNode=nodeGroup->begin(); itNode!=nodeGroup->end();itNode++)
    {
        try
        {
            NodeInternalForce(itNode->second, nodeForceLocal);
            if (nodeForceLocal.rows()!=rNodeForce.rows())
                throw MechanicsException(__PRETTY_FUNCTION__, "The number of displacement components is not equal for all members of the group.");
            rNodeForce+=nodeForceLocal;
        }
        catch(NuTo::MechanicsException & b)
        {
            b.AddMessage("[NuTo::StructureBase::NodeGroupInternalForce] Error getting gradient of internal potential.");
            throw;
        }
        catch(...)
        {
            throw MechanicsException(__PRETTY_FUNCTION__, "Error getting gradient of internal potential (unspecified exception).");
        }
    }
}

void NuTo::StructureBase::NodeInternalForce(const NodeBase* rNodePtr, Eigen::VectorXd& rNodeForce)
{
    try
    {
        std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>> elementOutputMap;
        elementOutputMap[Element::eOutput::INTERNAL_GRADIENT] = std::make_shared<ElementOutputBlockVectorDouble>(mDofStatus);
        elementOutputMap[Element::eOutput::GLOBAL_ROW_DOF] = std::make_shared<ElementOutputBlockVectorInt>(mDofStatus);

        std::vector<ElementBase*> elements;
        this->NodeGetElements(rNodePtr, elements);

        rNodeForce.resize(rNodePtr->GetNum(Node::eDof::DISPLACEMENTS));
        rNodeForce.setZero();

        for (auto element : elements)
        {
            element->Evaluate(elementOutputMap);
            const auto& internalGradient = elementOutputMap.at(Element::eOutput::INTERNAL_GRADIENT)->GetBlockFullVectorDouble()[Node::eDof::DISPLACEMENTS];
            const auto& globalRowDof = elementOutputMap.at(Element::eOutput::GLOBAL_ROW_DOF)->GetBlockFullVectorInt()[Node::eDof::DISPLACEMENTS];
            assert(internalGradient.rows() == globalRowDof.rows());

            for (int countDof=0; countDof< rNodePtr->GetNum(Node::eDof::DISPLACEMENTS); countDof++)
            {
                int theDof = rNodePtr->GetDof(Node::eDof::DISPLACEMENTS, countDof);
                for (int iDof=0; iDof < globalRowDof.rows(); iDof++)
                {
                    if (globalRowDof[iDof] == theDof)
                    {
                        rNodeForce(countDof)+=internalGradient(iDof);
                    }
                }
            }
        }


    }    catch(NuTo::MechanicsException & b)
    {
        b.AddMessage(std::string("[") + __PRETTY_FUNCTION__ + "] Error getting gradient of internal potential.");
        throw;
    }
    catch(...)
    {
        throw MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] Error getting gradient of internal potential (unspecified exception).");
    }
}

void NuTo::StructureBase::NodeGetElements(const int rNodeId, std::vector<int>& rElementNumbers)
{
    throw MechanicsException(__PRETTY_FUNCTION__, "Not available for this structure type.");
}

void NuTo::StructureBase::NodeGetElements(const NuTo::NodeBase* rNodePtr, std::vector<NuTo::ElementBase*>& rElements)
{
    throw MechanicsException(__PRETTY_FUNCTION__, "Not available for this structure type.");
}


int NuTo::StructureBase::NodeGetIdAtCoordinate(Eigen::VectorXd rCoordinates, double rRange)
{
    NuTo::Timer(__FUNCTION__, GetShowTime(), GetLogger());

    std::vector<std::pair<int,NodeBase*> > nodeVector;
    this->GetNodesTotal(nodeVector);

    double distance;

    int nodeId = -1;
    for (unsigned int countNode=0; countNode<nodeVector.size(); countNode++)
    {
        NodeBase* nodePtr(nodeVector[countNode].second);
        if (nodePtr->GetNum(Node::eDof::COORDINATES)<1)
            continue;

        distance = (nodePtr->Get(Node::eDof::COORDINATES)-rCoordinates).norm();

        if (distance<rRange)
        {
            if (nodeId==-1)
            {
                nodeId = nodeVector[countNode].first;
            }
            else
                throw MechanicsException(__PRETTY_FUNCTION__, "there is more than one node at that coordinate position.");
        }
    }
    if (nodeId==-1)
    {
        mLogger << "[NuTo::StructureBase::NodeGetIdAtCoordinate] no node could be found, return -1 as node id\n";
    }
    return nodeId;
}


#ifdef ENABLE_VISUALIZE
void NuTo::StructureBase::NodeTotalAddToVisualize(VisualizeUnstructuredGrid& rVisualize, const std::list<std::shared_ptr<NuTo::VisualizeComponent>>& rVisualizationList) const
{
    std::vector<const NodeBase*> nodeVec;
    this->GetNodesTotal(nodeVec);
    NodeVectorAddToVisualize(rVisualize,rVisualizationList,nodeVec);
}

void NuTo::StructureBase::NodeVectorAddToVisualize(VisualizeUnstructuredGrid& rVisualize, const std::list<std::shared_ptr<NuTo::VisualizeComponent>>& rVisualizationList, const std::vector<const NodeBase*>& rNodes) const
{
    for (unsigned int nodeCount = 0; nodeCount < rNodes.size(); nodeCount++)
    {
        rNodes[nodeCount]->Visualize(rVisualize, rVisualizationList);
    }
}
#endif //ENABLE_VISUALIZE
