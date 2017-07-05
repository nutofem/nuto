#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/loads/LoadNodeForces2D.h"

using namespace NuTo;

LoadNodeForces2D::LoadNodeForces2D(const NodeBase* rNode, const Eigen::MatrixXd& rDirection, double rValue)
    : LoadNode(rNode)
{
    if (rDirection.cols() != 1 || rDirection.rows() != 2)
        throw MechanicsException(__PRETTY_FUNCTION__,
                                 "Dimension of the direction matrix must be equal to the dimension of the structure.");

    memcpy(mDirection, rDirection.data(), 2 * sizeof(double));
    // normalize the direction
    double norm = sqrt(mDirection[0] * mDirection[0] + mDirection[1] * mDirection[1]);
    if (norm < 1e-14)
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Direction vector has zero length");
    }
    double invNorm = 1. / norm;
    mDirection[0] *= invNorm;
    mDirection[1] *= invNorm;
    mValue = rValue;
}


void LoadNodeForces2D::AddLoadToGlobalSubVectors(StructureOutputBlockVector& externalLoad) const
{
    assert(externalLoad.J[Node::eDof::DISPLACEMENTS].cols() == 1);
    assert(externalLoad.K[Node::eDof::DISPLACEMENTS].cols() == 1);
    for (int dofCount = 0; dofCount < 2; dofCount++)
    {
        int dof = mNode->GetDof(Node::eDof::DISPLACEMENTS, dofCount);
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

