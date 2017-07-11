#include "base/Exception.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/groups/Group.h"
#include "mechanics/loads/LoadNodeGroupForces1D.h"

using namespace NuTo;

LoadNodeGroupForces1D::LoadNodeGroupForces1D(const Group<NodeBase>* rGroup, double rDirection, double rValue)
    : LoadNodeGroup(rGroup)
{
    // set direction
    if (std::abs(rDirection) < 1e-14)
    {
        throw Exception(__PRETTY_FUNCTION__, "Direction vector has zero length.");
    }
    mDirection = rDirection / std::abs(rDirection);

    // set value
    mValue = rValue;
}


void LoadNodeGroupForces1D::AddLoadToGlobalSubVectors(StructureOutputBlockVector& externalLoad) const
{
    assert(externalLoad.J[Node::eDof::DISPLACEMENTS].cols() == 1);
    assert(externalLoad.K[Node::eDof::DISPLACEMENTS].cols() == 1);
    for (auto node : *mGroup)
    {
        int dof = node.second->GetDof(Node::eDof::DISPLACEMENTS, 0);
        assert(dof >= 0);
        if (dof < externalLoad.J[Node::eDof::DISPLACEMENTS].rows())
        {
            externalLoad.J[Node::eDof::DISPLACEMENTS](dof, 0) += mDirection * mValue;
        }
        else
        {
            dof -= externalLoad.J[Node::eDof::DISPLACEMENTS].rows();
            assert(dof < externalLoad.K[Node::eDof::DISPLACEMENTS].rows());
            externalLoad.K[Node::eDof::DISPLACEMENTS](dof, 0) += mDirection * mValue;
        }
    }
}
