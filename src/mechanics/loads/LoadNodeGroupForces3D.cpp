#include "base/Exception.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/groups/Group.h"
#include "mechanics/loads/LoadNodeGroupForces3D.h"

using namespace NuTo;

LoadNodeGroupForces3D::LoadNodeGroupForces3D(const Group<NodeBase>* rGroup, const Eigen::VectorXd& rDirection,
                                             double rValue)
    : LoadNodeGroup(rGroup)
{
    if (rDirection.rows() != 3)
        throw Exception(__PRETTY_FUNCTION__,
                                 "Dimension of the direction matrix must be equal to the dimension of the structure.");

    memcpy(mDirection, rDirection.data(), 3 * sizeof(double));
    // normalize the direction
    double norm = sqrt(mDirection[0] * mDirection[0] + mDirection[1] * mDirection[1] + mDirection[2] * mDirection[2]);
    if (norm < 1e-14)
    {
        throw Exception(__PRETTY_FUNCTION__, "Direction vector has zero length");
    }
    double invNorm = 1. / norm;
    mDirection[0] *= invNorm;
    mDirection[1] *= invNorm;
    mDirection[2] *= invNorm;
    mValue = rValue;
}


void LoadNodeGroupForces3D::AddLoadToGlobalSubVectors(StructureOutputBlockVector& externalLoad) const
{
    assert(externalLoad.J[Node::eDof::DISPLACEMENTS].cols() == 1);
    assert(externalLoad.K[Node::eDof::DISPLACEMENTS].cols() == 1);
    for (auto node : *mGroup)
    {
        for (int dofCount = 0; dofCount < 3; dofCount++)
        {
            int dof = node.second->GetDof(Node::eDof::DISPLACEMENTS, dofCount);
            assert(dof >= 0);
            if (dof < externalLoad.J[Node::eDof::DISPLACEMENTS].rows())
            {
                externalLoad.J[Node::eDof::DISPLACEMENTS](dof, 0) += this->mValue * mDirection[dofCount];
            }
            else
            {
                dof -= externalLoad.J[Node::eDof::DISPLACEMENTS].rows();
                assert(dof < externalLoad.K[Node::eDof::DISPLACEMENTS].rows());
                externalLoad.K[Node::eDof::DISPLACEMENTS](dof, 0) += this->mValue * mDirection[dofCount];
            }
        }
    }
}
