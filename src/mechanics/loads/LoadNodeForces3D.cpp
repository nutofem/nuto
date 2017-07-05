#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/loads/LoadNodeForces3D.h"

using namespace NuTo;

LoadNodeForces3D::LoadNodeForces3D(const NodeBase* rNode, const Eigen::MatrixXd& rDirection, double rValue)
    : LoadNode(rNode)
{
    if (rDirection.cols() != 1 || rDirection.rows() != 3)
        throw MechanicsException(__PRETTY_FUNCTION__,
                                 "Dimension of the direction matrix must be equal to the dimension of the structure.");

    memcpy(mDirection, rDirection.data(), 3 * sizeof(double));
    // normalize the direction
    double norm = sqrt(mDirection[0] * mDirection[0] + mDirection[1] * mDirection[1] + mDirection[2] * mDirection[2]);
    if (norm < 1e-14)
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Direction vector has zero length");
    }
    double invNorm = 1. / norm;
    mDirection[0] *= invNorm;
    mDirection[1] *= invNorm;
    mDirection[2] *= invNorm;
    mValue = rValue;
}


void LoadNodeForces3D::AddLoadToGlobalSubVectors(StructureOutputBlockVector& externalLoad) const
{
    assert(externalLoad.J[Node::eDof::DISPLACEMENTS].cols() == 1);
    assert(externalLoad.K[Node::eDof::DISPLACEMENTS].cols() == 1);
    for (int dofCount = 0; dofCount < 3; dofCount++)
    {
        int dof = mNode->GetDof(Node::eDof::DISPLACEMENTS, dofCount);
        assert(dof >= 0);
        if (dof < externalLoad.J[Node::eDof::DISPLACEMENTS].rows())
        {
            externalLoad.J[Node::eDof::DISPLACEMENTS](dof, 0) += mValue * mDirection[dofCount];
        }
        else
        {
            dof -= externalLoad.J[Node::eDof::DISPLACEMENTS].rows();
            assert(dof < externalLoad.K[Node::eDof::DISPLACEMENTS].rows());
            externalLoad.K[Node::eDof::DISPLACEMENTS](dof, 0) += mValue * mDirection[dofCount];
        }
    }
}
