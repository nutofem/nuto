#include "base/Timer.h"

#include "math/EigenCompanion.h"

#include "mechanics/elements/ElementOutputBlockVectorDouble.h"
#include "mechanics/elements/ElementOutputBlockVectorInt.h"

#include "mechanics/structures/StructureBase.h"
#include "mechanics/structures/StructureOutputBlockVector.h"
#include "mechanics/elements/ElementBase.h"
#include "mechanics/elements/ElementEnum.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/groups/Group.h"
#include "mechanics/groups/GroupEnum.h"

#ifdef ENABLE_VISUALIZE
#include "visualize/UnstructuredGrid.h"
#include "visualize/ComponentName.h"
#endif

void NuTo::StructureBase::NodeSetDisplacements(int rNode, const Eigen::VectorXd& rDisplacements)
{
    NuTo::Timer(__FUNCTION__, GetShowTime(), GetLogger());
    NodeBase* nodePtr = NodeGetNodePtr(rNode);
    this->mUpdateTmpStaticDataRequired = true;

    if (rDisplacements.cols() != 1)
        throw Exception(__PRETTY_FUNCTION__, "Displacement matrix has to have a single column.");

    if (rDisplacements.rows() <= 0 or rDisplacements.rows() > 3)
        throw Exception(__PRETTY_FUNCTION__, "The number of displacement components is either 1, 2 or 3.");

    nodePtr->Set(Node::eDof::DISPLACEMENTS, rDisplacements);
}


NuTo::StructureOutputBlockVector NuTo::StructureBase::NodeExtractDofValues() const
{
    return NodeExtractDofValues(0);
}

void NuTo::StructureBase::NodeSetDisplacements(int rNode, int rTimeDerivative, const Eigen::VectorXd& rDisplacements)
{
    NuTo::Timer(__FUNCTION__, GetShowTime(), GetLogger());
    NodeBase* nodePtr = NodeGetNodePtr(rNode);
    this->mUpdateTmpStaticDataRequired = true;

    if (rDisplacements.cols() != 1)
        throw Exception(__PRETTY_FUNCTION__, "Displacement matrix has to have a single column.");
    if (nodePtr->GetNumTimeDerivatives(Node::eDof::DISPLACEMENTS) < rTimeDerivative)
        throw Exception(__PRETTY_FUNCTION__,
                                 "number of time derivatives stored at node is less than the required value.");
    if (rDisplacements.rows() <= 0 or rDisplacements.rows() > 3)
        throw Exception(__PRETTY_FUNCTION__, "The number of displacement components is either 1, 2 or 3.");

    nodePtr->Set(Node::eDof::DISPLACEMENTS, rTimeDerivative, rDisplacements);
}

void NuTo::StructureBase::NodeGroupSetDisplacements(int rGroupIdent, const Eigen::VectorXd& rDisplacements)
{
    NuTo::Timer(__FUNCTION__, GetShowTime(), GetLogger());
    this->mUpdateTmpStaticDataRequired = true;
    if (rDisplacements.cols() != 1)
        throw Exception(__PRETTY_FUNCTION__, "Displacement matrix has to have a single column.");

    boost::ptr_map<int, GroupBase>::iterator itGroup = mGroupMap.find(rGroupIdent);
    if (itGroup == mGroupMap.end())
        throw Exception(__PRETTY_FUNCTION__, "Group with the given identifier does not exist.");
    if (itGroup->second->GetType() != NuTo::eGroupId::Nodes)
        throw Exception(__PRETTY_FUNCTION__, "Group is not a node group.");
    Group<NodeBase>* nodeGroup = dynamic_cast<Group<NodeBase>*>(itGroup->second);
    assert(nodeGroup != 0);

    for (auto& node : *nodeGroup)
    {
        if (rDisplacements.rows() <= 0 or rDisplacements.rows() > 3)
            throw Exception(__PRETTY_FUNCTION__, "The number of displacement components is either 1, 2 or 3.");

        node.second->Set(Node::eDof::DISPLACEMENTS, rDisplacements);
    }
}

void NuTo::StructureBase::NodeGroupSetDisplacements(int rGroupIdent, int rTimeDerivative,
                                                    const Eigen::VectorXd& rDisplacements)
{
    NuTo::Timer(__FUNCTION__, GetShowTime(), GetLogger());
    this->mUpdateTmpStaticDataRequired = true;
    if (rDisplacements.cols() != 1)
        throw Exception(__PRETTY_FUNCTION__, "Displacement matrix has to have a single column.");

    boost::ptr_map<int, GroupBase>::iterator itGroup = mGroupMap.find(rGroupIdent);
    if (itGroup == mGroupMap.end())
        throw Exception(__PRETTY_FUNCTION__, "Group with the given identifier does not exist.");
    if (itGroup->second->GetType() != NuTo::eGroupId::Nodes)
        throw Exception(__PRETTY_FUNCTION__, "Group is not a node group.");
    Group<NodeBase>* nodeGroup = dynamic_cast<Group<NodeBase>*>(itGroup->second);
    assert(nodeGroup != 0);

    for (auto& node : *nodeGroup)
    {
        if (node.second->GetNumTimeDerivatives(Node::eDof::DISPLACEMENTS) < rTimeDerivative)
            throw Exception(__PRETTY_FUNCTION__, "does not have a sufficient number of time derivatives.");

        if (rDisplacements.rows() <= 0 or rDisplacements.rows() > 3)
            throw Exception(__PRETTY_FUNCTION__, "The number of displacement components is either 1, 2 or 3.");

        node.second->Set(Node::eDof::DISPLACEMENTS, rTimeDerivative, rDisplacements);
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
    if (nodePtr->GetNumTimeDerivatives(Node::eDof::TEMPERATURE) < rTimeDerivative)
        throw Exception(__PRETTY_FUNCTION__,
                                 "Number of time derivatives stored at node is less than the required value.");
    nodePtr->Set(Node::eDof::TEMPERATURE, rTimeDerivative, rTemperature);
}

void NuTo::StructureBase::NodeGroupGetMembers(int rGroupId, std::vector<int>& rMembers)
{
    NuTo::Timer(__FUNCTION__, GetShowTime(), GetLogger());

    boost::ptr_map<int, GroupBase>::iterator itGroup = mGroupMap.find(rGroupId);
    if (itGroup == mGroupMap.end())
        throw Exception(__PRETTY_FUNCTION__, "Group with the given identifier does not exist.");
    if (itGroup->second->GetType() != NuTo::eGroupId::Nodes)
        throw Exception(__PRETTY_FUNCTION__, "Group is not an node group.");
    Group<NodeBase>* nodeGroup = itGroup->second->AsGroupNode();
    assert(nodeGroup != 0);

    rMembers.resize(nodeGroup->GetNumMembers());
    int countNode(0);
    for (Group<NodeBase>::const_iterator itNode = nodeGroup->begin(); itNode != nodeGroup->end(); ++itNode, countNode++)
    {
        rMembers[countNode] = itNode->first;
    }
}


void NuTo::StructureBase::NodeGetDisplacements(int rNode, Eigen::VectorXd& rDisplacements) const
{
    this->NodeGetDisplacements(rNode, 0, rDisplacements);
}


void NuTo::StructureBase::NodeGetDisplacements(int rNode, int rTimeDerivative, Eigen::VectorXd& rDisplacements) const
{
    NuTo::Timer(__FUNCTION__, GetShowTime(), GetLogger());

    const NodeBase* nodePtr = NodeGetNodePtr(rNode);

    if (nodePtr->GetNum(Node::eDof::DISPLACEMENTS) == 0)
        throw Exception(__PRETTY_FUNCTION__, "Node has no displacements.");

    rDisplacements = nodePtr->Get(Node::eDof::DISPLACEMENTS, rTimeDerivative);
}


std::vector<int> NuTo::StructureBase::NodeGetDofIds(const int rNodeId, NuTo::Node::eDof rDof) const
{

    const NodeBase* nodePtr = NodeGetNodePtr(rNodeId);

    int numDofIds = nodePtr->GetNum(rDof);

    if (numDofIds == 0)
        throw Exception(__PRETTY_FUNCTION__, "Node does not have the requested dof.");

    std::vector<int> dofIds(numDofIds);
    for (int i = 0; i < numDofIds; ++i)
        dofIds[i] = nodePtr->GetDof(rDof, i);

    return dofIds;
}


void NuTo::StructureBase::NodeGroupGetDisplacements(int rGroupIdent, Eigen::MatrixXd& rDisplacements)
{
    NuTo::Timer(__FUNCTION__, GetShowTime(), GetLogger());

    boost::ptr_map<int, GroupBase>::iterator itGroup = mGroupMap.find(rGroupIdent);
    if (itGroup == mGroupMap.end())
        throw Exception(__PRETTY_FUNCTION__, "Group with the given identifier does not exist.");
    if (itGroup->second->GetType() != NuTo::eGroupId::Nodes)
        throw Exception(__PRETTY_FUNCTION__, "Group is not a node group.");
    Group<NodeBase>* nodeGroup = itGroup->second->AsGroupNode();
    assert(nodeGroup != 0);

    // all nodes have to have the same dimension
    if (nodeGroup->GetNumMembers() < 1)
        throw Exception(__PRETTY_FUNCTION__, "Group has no members.");

    int numDisp = nodeGroup->begin()->second->GetNum(Node::eDof::DISPLACEMENTS);
    // resize the matrix
    rDisplacements.resize(nodeGroup->GetNumMembers(), numDisp);

    int theNode(0);
    for (Group<NodeBase>::iterator itNode = nodeGroup->begin(); itNode != nodeGroup->end(); ++itNode, theNode++)
    {
        if (numDisp != 1 and numDisp != 2 and numDisp != 3)
            throw Exception(__PRETTY_FUNCTION__,
                                     "The number of displacement components is either 1, 2 or 3.");

        rDisplacements.row(theNode) = itNode->second->Get(Node::eDof::DISPLACEMENTS).transpose();
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
        throw Exception(__PRETTY_FUNCTION__, "Node doesn't have a temperature.");
    return nodePtr->Get(Node::eDof::TEMPERATURE, rTimeDerivative)[0];
}

void NuTo::StructureBase::NodeGetCoordinates(int rNode, Eigen::VectorXd& rCoordinates) const
{
    NuTo::Timer(__FUNCTION__, GetShowTime(), GetLogger());

    const NodeBase* nodePtr = NodeGetNodePtr(rNode);

    if (nodePtr->GetNum(Node::eDof::COORDINATES) == 0)
        throw Exception(__PRETTY_FUNCTION__, "Node has no coordinates.");

    rCoordinates = nodePtr->Get(Node::eDof::COORDINATES);
}

void NuTo::StructureBase::NodeGroupGetCoordinates(int rGroupIdent, Eigen::MatrixXd& rCoordinates)
{
    NuTo::Timer(__FUNCTION__, GetShowTime(), GetLogger());

    boost::ptr_map<int, GroupBase>::iterator itGroup = mGroupMap.find(rGroupIdent);
    if (itGroup == mGroupMap.end())
        throw Exception(__PRETTY_FUNCTION__, "Group with the given identifier does not exist.");
    if (itGroup->second->GetType() != NuTo::eGroupId::Nodes)
        throw Exception(__PRETTY_FUNCTION__, "Group is not a node group.");
    Group<NodeBase>* nodeGroup = itGroup->second->AsGroupNode();
    assert(nodeGroup != 0);

    // all nodes have to have the same dimension
    if (nodeGroup->GetNumMembers() < 1)
        throw Exception(__PRETTY_FUNCTION__, "Group has no members.");

    int numCoords = nodeGroup->begin()->second->GetNum(Node::eDof::COORDINATES);
    // resize the matrix
    rCoordinates.resize(nodeGroup->GetNumMembers(), numCoords);
    int theNode(0);
    for (Group<NodeBase>::iterator itNode = nodeGroup->begin(); itNode != nodeGroup->end(); ++itNode, theNode++)
    {
        if (numCoords != 1 and numCoords != 2 and numCoords != 3)
            throw Exception(__PRETTY_FUNCTION__,
                                     "The number of coordinates components is either 1, 2 or 3.");

        rCoordinates.row(theNode) = itNode->second->Get(Node::eDof::COORDINATES).transpose();
    }
}


void NuTo::StructureBase::NodeInternalForce(int rId, Eigen::VectorXd& rNodeForce)
{
    NuTo::Timer(__FUNCTION__, GetShowTime(), GetLogger());

    const NodeBase* nodePtr = NodeGetNodePtr(rId);
    NodeInternalForce(nodePtr, rNodeForce);
}

void NuTo::StructureBase::NodeGroupInternalForce(int rGroupIdent, Eigen::VectorXd& rNodeForce)
{
    NuTo::Timer(__FUNCTION__, GetShowTime(), GetLogger());

    boost::ptr_map<int, GroupBase>::const_iterator itGroup = mGroupMap.find(rGroupIdent);
    if (itGroup == mGroupMap.end())
        throw Exception(__PRETTY_FUNCTION__, "Group with the given identifier does not exist.");
    if (itGroup->second->GetType() != NuTo::eGroupId::Nodes)
        throw Exception(__PRETTY_FUNCTION__, "Group is not a node group.");
    const Group<NodeBase>* nodeGroup = dynamic_cast<const Group<NodeBase>*>(itGroup->second);
    assert(nodeGroup != 0);

    Eigen::VectorXd nodeForceLocal;

    if (nodeGroup->GetNumMembers() == 0)
        throw Exception(__PRETTY_FUNCTION__, "Node group is empty.");
    rNodeForce.resize(nodeGroup->begin()->second->GetNum(Node::eDof::DISPLACEMENTS));
    rNodeForce.setZero();

    for (auto node : *nodeGroup)
    {
        NodeInternalForce(node.second, nodeForceLocal);
        if (nodeForceLocal.rows() != rNodeForce.rows())
            throw Exception(
                    __PRETTY_FUNCTION__,
                    "The number of displacement components is not equal for all members of the group.");
        rNodeForce += nodeForceLocal;
    }
}

void NuTo::StructureBase::NodeInternalForce(const NodeBase* rNodePtr, Eigen::VectorXd& rNodeForce)
{
    std::map<ElementEnum::eOutput, std::shared_ptr<ElementOutputBase>> elementOutputMap;
    elementOutputMap[ElementEnum::eOutput::INTERNAL_GRADIENT] =
            std::make_shared<ElementOutputBlockVectorDouble>(GetDofStatus());
    elementOutputMap[ElementEnum::eOutput::GLOBAL_ROW_DOF] =
            std::make_shared<ElementOutputBlockVectorInt>(GetDofStatus());

    std::vector<ElementBase*> elements;
    this->NodeGetElements(rNodePtr, elements);

    rNodeForce.resize(rNodePtr->GetNum(Node::eDof::DISPLACEMENTS));
    rNodeForce.setZero();

    for (auto element : elements)
    {
        element->Evaluate(elementOutputMap);
        const auto& internalGradient = elementOutputMap.at(ElementEnum::eOutput::INTERNAL_GRADIENT)
                                               ->GetBlockFullVectorDouble()[Node::eDof::DISPLACEMENTS];
        const auto& globalRowDof = elementOutputMap.at(ElementEnum::eOutput::GLOBAL_ROW_DOF)
                                           ->GetBlockFullVectorInt()[Node::eDof::DISPLACEMENTS];
        assert(internalGradient.rows() == globalRowDof.rows());

        for (int countDof = 0; countDof < rNodePtr->GetNum(Node::eDof::DISPLACEMENTS); countDof++)
        {
            int theDof = rNodePtr->GetDof(Node::eDof::DISPLACEMENTS, countDof);
            for (int iDof = 0; iDof < globalRowDof.rows(); iDof++)
            {
                if (globalRowDof[iDof] == theDof)
                {
                    rNodeForce(countDof) += internalGradient(iDof);
                }
            }
        }
    }
}

void NuTo::StructureBase::NodeGetElements(const int, std::vector<int>&)
{
    throw Exception(__PRETTY_FUNCTION__, "Not available for this structure type.");
}

void NuTo::StructureBase::NodeGetElements(const NuTo::NodeBase*, std::vector<NuTo::ElementBase*>&)
{
    throw Exception(__PRETTY_FUNCTION__, "Not available for this structure type.");
}

NuTo::NodeBase& NuTo::StructureBase::NodeGetAtCoordinate(double coordinate, double tolerance)
{
    if (GetDimension() != 1)
        throw Exception(__PRETTY_FUNCTION__, "Only valid for 1D structure!");
    return NodeGetAtCoordinate(Eigen::Matrix<double, 1, 1>::Constant(coordinate), tolerance);
}

NuTo::NodeBase& NuTo::StructureBase::NodeGetAtCoordinate(Eigen::VectorXd coordinate, double tolerance)
{
    NuTo::Timer(__FUNCTION__, GetShowTime(), GetLogger());

    std::vector<NodeBase*> nodeVector;
    this->GetNodesTotal(nodeVector);

    double toleranceSquared = tolerance * tolerance;

    for (auto* node : nodeVector)
    {
        if (node->GetNum(Node::eDof::COORDINATES) < 1)
            continue;
        if ((node->Get(Node::eDof::COORDINATES) - coordinate).squaredNorm() < toleranceSquared)
            return *node;
    }
    std::stringstream coordStream;
    coordStream << '(' << coordinate.transpose() << ')';
    throw Exception(__PRETTY_FUNCTION__, "There is no node at " + coordStream.str() + " within tolerance " +
                                                          std::to_string(tolerance));
}

int NuTo::StructureBase::NodeGetIdAtCoordinate(Eigen::VectorXd rCoordinates, double rRange)
{
    NuTo::Timer(__FUNCTION__, GetShowTime(), GetLogger());

    std::vector<std::pair<int, NodeBase*>> nodeVector;
    this->GetNodesTotal(nodeVector);

    double distance;

    int nodeId = -1;
    for (auto& node : nodeVector)
    {
        NodeBase* nodePtr(node.second);
        if (nodePtr->GetNum(Node::eDof::COORDINATES) < 1)
            continue;

        distance = (nodePtr->Get(Node::eDof::COORDINATES) - rCoordinates).norm();

        if (distance < rRange)
        {
            if (nodeId == -1)
            {
                nodeId = node.first;
            }
            else
                throw Exception(__PRETTY_FUNCTION__,
                                         "there is more than one node at that coordinate position.");
        }
    }
    if (nodeId == -1)
    {
        mLogger << "[NuTo::StructureBase::NodeGetIdAtCoordinate] no node could be found, return -1 as node id\n";
    }
    return nodeId;
}


#ifdef ENABLE_VISUALIZE
void NuTo::StructureBase::NodeTotalAddToVisualize(Visualize::UnstructuredGrid& visualizer,
                                                  const std::vector<eVisualizeWhat>& visualizeComponents) const
{
    std::vector<const NodeBase*> nodeVec;
    this->GetNodesTotal(nodeVec);
    NodeVectorAddToVisualize(visualizer, visualizeComponents, nodeVec);
}

void NuTo::StructureBase::NodeVectorAddToVisualize(Visualize::UnstructuredGrid& visualizer,
                                                   const std::vector<eVisualizeWhat>& visualizeComponents,
                                                   const std::vector<const NodeBase*>& nodes) const
{
    using Node::eDof;

    auto NodeData3D = [=](const NodeBase* node, eDof dof, int timeDerivative) {
        auto data = node->Get(dof, timeDerivative);
        return EigenCompanion::To3D(data);
    };

    for (auto node : nodes)
    {
        /*
         * TTitscher: vtk requires each point to have all data fields. Imagine:
         * Node 0 : Displ, Temp
         * Node 1 : Displ
         * Node 2 : Displ, Temp, Rotation
         *
         * Each node is a point. "Displ", "Temp" and "Rotation" are data fields.
         * What value should be assigned to "Rotation" of Node 1? --> Problem.
         *
         * Near future solution: Independent interpolation, nodes have only one dof
         * --> NodesDispl.vtu, NodesTemp.vtu, NodesRotation.vtu.
         *  TODO for future self.
         *
         *  For now, only displacements (+derivatives) are visualized (as before...).
         */

        if (node->GetNum(Node::eDof::DISPLACEMENTS) == 0)
            continue;
        auto coordinates = NuTo::EigenCompanion::To3D(node->Get(eDof::COORDINATES));
        const int pointId = visualizer.AddPoint(coordinates);


        // store data
        for (auto component : visualizeComponents)
        {
            auto name = GetComponentName(component);
            switch (component)
            {
            case eVisualizeWhat::DISPLACEMENTS:
            {
                visualizer.SetPointData(pointId, name, NodeData3D(node, Node::eDof::DISPLACEMENTS, 0));
            }
            break;
            case eVisualizeWhat::VELOCITY:
            {
                if (node->GetNumTimeDerivatives(Node::eDof::DISPLACEMENTS) < 1)
                    break;
                visualizer.SetPointData(pointId, name, NodeData3D(node, Node::eDof::DISPLACEMENTS, 1));
            }
            break;
            case eVisualizeWhat::ACCELERATION:
            {
                if (node->GetNumTimeDerivatives(Node::eDof::DISPLACEMENTS) < 2)
                    break;
                visualizer.SetPointData(pointId, name, NodeData3D(node, Node::eDof::DISPLACEMENTS, 2));
            }
            break;
            default:
                break;
            }
        }
    }
}
#endif // ENABLE_VISUALIZE
