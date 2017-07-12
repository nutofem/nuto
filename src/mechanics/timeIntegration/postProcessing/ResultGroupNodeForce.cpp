#include "mechanics/timeIntegration/postProcessing/ResultGroupNodeForce.h"

#include "mechanics/structures/StructureBase.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/groups/Group.h"

using namespace NuTo;

ResultGroupNodeForce::ResultGroupNodeForce(const std::string& rIdent, int rGroupNodeId)
    : ResultGroupNodeDof(rIdent, rGroupNodeId)
{
}


int ResultGroupNodeForce::GetNumData(const StructureBase& rStructure) const
{
    const Group<NodeBase>& groupNode = *rStructure.GroupGetGroupPtr(mGroupNodeId)->AsGroupNode();

    // all nodes have to have the same dimension (number of displacement components)
    if (groupNode.GetNumMembers() < 1)
        throw MechanicsException("[ResultGroupNodeForce::GetNumData] Group has no members.");

    return groupNode.begin()->second->GetNum(Node::eDof::DISPLACEMENTS);
}

Eigen::VectorXd ResultGroupNodeForce::CalculateValues(const StructureBase& rStructure,
                                                      const StructureOutputBlockVector& residual) const
{
    const Group<NodeBase>& nodeGroup = *GetGroupNodePtr(rStructure);
    const int dim = rStructure.GetDimension();
    Eigen::VectorXd result = Eigen::VectorXd::Zero(dim);

    for (auto& itNode : nodeGroup)
    {
        assert(itNode.second->GetNum(Node::eDof::DISPLACEMENTS) == dim);
        for (int iDim = 0; iDim < dim; iDim++)
        {
            int theDof = itNode.second->GetDof(Node::eDof::DISPLACEMENTS, iDim);
            if (theDof < rStructure.GetNumActiveDofs(Node::eDof::DISPLACEMENTS))
            {
                result[iDim] -= residual.J[Node::eDof::DISPLACEMENTS](theDof); // Defined as minus. [R = F_ext - F_int]
            }
            else
            {
                result[iDim] -= residual.K[Node::eDof::DISPLACEMENTS](theDof - rStructure.GetNumActiveDofs(Node::eDof::DISPLACEMENTS));
            }
        }
    }
    return result;
}
