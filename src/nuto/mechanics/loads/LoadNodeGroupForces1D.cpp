// $Id$
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/loads/LoadNodeGroupForces1D.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"

//! @brief constructor
NuTo::LoadNodeGroupForces1D::LoadNodeGroupForces1D(int rLoadCase, const Group<NodeBase>* rGroup, double rDirection, double rValue) :
        LoadNodeGroup(rLoadCase,rGroup)
{
    // set direction
    if (fabs(rDirection) < 1e-14)
    {
        throw MechanicsException("[NuTo::LoadNodeGroupForces1D::LoadNodeGroupForces1D] direction vector has zero length.");
    }
    this->mDirection = rDirection/fabs(rDirection);

    // set value
    this->mValue = rValue;
}

// adds the load to global sub-vectors
void NuTo::LoadNodeGroupForces1D::AddLoadToGlobalSubVectors(int rLoadCase, NuTo::FullVector<double,Eigen::Dynamic>& rActiceDofsLoadVector, NuTo::FullVector<double,Eigen::Dynamic>& rDependentDofsLoadVector)const
{
    if (rLoadCase!=mLoadCase)
    	return;
    assert(rActiceDofsLoadVector.GetNumColumns()==1);
    assert(rDependentDofsLoadVector.GetNumColumns()==1);
    for (Group<NodeBase>::const_iterator itNode=this->mGroup->begin(); itNode!=this->mGroup->end(); itNode++)
    {
        try
        {
            int dof = itNode->second->GetDofDisplacement(0);
            assert(dof >= 0);
            if (dof < rActiceDofsLoadVector.GetNumRows())
            {
                rActiceDofsLoadVector(dof,0) += this->mDirection * this->mValue;
            }
            else
            {
                dof -= rActiceDofsLoadVector.GetNumRows();
                assert(dof < rDependentDofsLoadVector.GetNumRows());
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
