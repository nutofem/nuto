// $Id$
#include "mechanics/MechanicsException.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/groups/Group.h"
#include "mechanics/loads/LoadNodeGroupForces1D.h"

#include "math/SparseMatrixCSRGeneral.h"

//! @brief constructor
NuTo::LoadNodeGroupForces1D::LoadNodeGroupForces1D(int rLoadCase, const Group<NodeBase>* rGroup, double rDirection, double rValue) :
        LoadNodeGroup(rLoadCase,rGroup)
{
    // set direction
    if (std::abs(rDirection) < 1e-14)
    {
        throw MechanicsException("[NuTo::LoadNodeGroupForces1D::LoadNodeGroupForces1D] direction vector has zero length.");
    }
    this->mDirection = rDirection/std::abs(rDirection);

    // set value
    this->mValue = rValue;
}

// adds the load to global sub-vectors
void NuTo::LoadNodeGroupForces1D::AddLoadToGlobalSubVectors(int rLoadCase, Eigen::VectorXd& rActiceDofsLoadVector, Eigen::VectorXd& rDependentDofsLoadVector)const
{
    if (rLoadCase!=mLoadCase)
    	return;
    assert(rActiceDofsLoadVector.cols()==1);
    assert(rDependentDofsLoadVector.cols()==1);
    for (Group<NodeBase>::const_iterator itNode=this->mGroup->begin(); itNode!=this->mGroup->end(); itNode++)
    {
        try
        {
            int dof = itNode->second->GetDof(Node::eDof::DISPLACEMENTS, 0);
            assert(dof >= 0);
            if (dof < rActiceDofsLoadVector.rows())
            {
                rActiceDofsLoadVector(dof,0) += this->mDirection * this->mValue;
            }
            else
            {
                dof -= rActiceDofsLoadVector.rows();
                assert(dof < rDependentDofsLoadVector.rows());
                rDependentDofsLoadVector(dof,0) += this->mDirection * this->mValue;
            }
        }
        catch (std::bad_cast & b)
        {
            throw MechanicsException("[NuTo::LoadNodeForces1D::AddLoad] Node has no displacements or its dimension is not equivalent to the 1.");
        }
        catch (...)
        {
            throw MechanicsException("[NuTo::LoadNodeForces1D::AddLoad] Error getting displacements of node (unspecified exception).");
        }
    }
}
