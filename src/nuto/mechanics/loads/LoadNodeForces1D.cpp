// $Id$
#include <assert.h>
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/nodes/NodeEnum.h"
#include "nuto/mechanics/loads/LoadNodeForces1D.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"

// constructor
NuTo::LoadNodeForces1D::LoadNodeForces1D(int rLoadCase, const NodeBase* rNode, double rDirection, double rValue) :
        LoadNode(rLoadCase,rNode)
{
    // set direction
    if (std::abs(rDirection) < 1e-14)
    {
        throw MechanicsException("[NuTo::LoadNodeForces1D::LoadNodeForces1D] direction vector has zero length.");
    }
    this->mDirection = rDirection/std::abs(rDirection);

    // set value
    this->mValue = rValue;
}

// adds the load to global sub-vectors
void NuTo::LoadNodeForces1D::AddLoadToGlobalSubVectors(int rLoadCase, NuTo::FullVector<double,Eigen::Dynamic>& rActiveDofsLoadVector, NuTo::FullVector<double,Eigen::Dynamic>& rDependentDofsLoadVector)const
{
    if (rLoadCase!=mLoadCase)
    	return;
    assert(rActiveDofsLoadVector.GetNumColumns()==1);
    assert(rDependentDofsLoadVector.GetNumColumns()==1);
    try
    {
        int dof = mNode->GetDof(Node::eDof::DISPLACEMENTS, 0);
        assert(dof >= 0);
        if (dof < rActiveDofsLoadVector.GetNumRows())
        {
            rActiveDofsLoadVector(dof,0) += this->mDirection * this->mValue;
        }
        else
        {
            dof -= rActiveDofsLoadVector.GetNumRows();
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

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::LoadNodeForces1D)
BOOST_CLASS_TRACKING(NuTo::LoadNodeForces1D, track_always)
#endif