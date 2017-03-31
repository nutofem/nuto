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


void LoadNodeForces1D::AddLoadToGlobalSubVectors(Eigen::VectorXd& rActiveDofsLoadVector,
                                                 Eigen::VectorXd& rDependentDofsLoadVector) const
{
    assert(rActiveDofsLoadVector.cols() == 1);
    assert(rDependentDofsLoadVector.cols() == 1);

    int dof = mNode->GetDof(Node::eDof::DISPLACEMENTS, 0);
    assert(dof >= 0);
    if (dof < rActiveDofsLoadVector.rows())
    {
        rActiveDofsLoadVector(dof, 0) += this->mDirection * this->mValue;
    }
    else
    {
        dof -= rActiveDofsLoadVector.rows();
        assert(dof < rDependentDofsLoadVector.rows());
        rDependentDofsLoadVector(dof, 0) += this->mDirection * this->mValue;
    }
}

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_IMPLEMENT(LoadNodeForces1D)
BOOST_CLASS_TRACKING(LoadNodeForces1D, track_always)
#endif
