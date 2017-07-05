#include <cassert>
#include "mechanics/MechanicsException.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/loads/LoadNodeForces1D.h"

using namespace NuTo;

LoadNodeForces1D::LoadNodeForces1D(const NodeBase* rNode, double rDirection, double rValue)
    : LoadNode(rNode)
{
    // set direction
    if (std::abs(rDirection) < 1e-14)
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Direction vector has zero length.");
    }
    mDirection = rDirection / std::abs(rDirection);

    // set value
    mValue = rValue;
}


void LoadNodeForces1D::AddLoadToGlobalSubVectors(StructureOutputBlockVector& externalLoad) const
{
    assert(externalLoad.J[Node::eDof::DISPLACEMENTS].cols() == 1);
    assert(externalLoad.K[Node::eDof::DISPLACEMENTS].cols() == 1);

    int dof = mNode->GetDof(Node::eDof::DISPLACEMENTS, 0);
    assert(dof >= 0);
    if (dof < externalLoad.J[Node::eDof::DISPLACEMENTS].rows())
    {
        externalLoad.J[Node::eDof::DISPLACEMENTS](dof, 0) += this->mDirection * this->mValue;
    }
    else
    {
        dof -= externalLoad.J[Node::eDof::DISPLACEMENTS].rows();
        assert(dof < externalLoad.K[Node::eDof::DISPLACEMENTS].rows());
        externalLoad.K[Node::eDof::DISPLACEMENTS](dof, 0) += this->mDirection * this->mValue;
    }
}

