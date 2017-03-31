#include "mechanics/MechanicsException.h"
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
        throw MechanicsException(__PRETTY_FUNCTION__, "Direction vector has zero length.");
    }
    mDirection = rDirection / std::abs(rDirection);

    // set value
    mValue = rValue;
}


void LoadNodeGroupForces1D::AddLoadToGlobalSubVectors(Eigen::VectorXd& rActiceDofsLoadVector,
                                                      Eigen::VectorXd& rDependentDofsLoadVector) const
{
    assert(rActiceDofsLoadVector.cols() == 1);
    assert(rDependentDofsLoadVector.cols() == 1);
    for (auto node : *mGroup)
    {
        int dof = node.second->GetDof(Node::eDof::DISPLACEMENTS, 0);
        assert(dof >= 0);
        if (dof < rActiceDofsLoadVector.rows())
        {
            rActiceDofsLoadVector(dof, 0) += mDirection * mValue;
        }
        else
        {
            dof -= rActiceDofsLoadVector.rows();
            assert(dof < rDependentDofsLoadVector.rows());
            rDependentDofsLoadVector(dof, 0) += mDirection * mValue;
        }
    }
}
