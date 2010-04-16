// $Id$
#include <assert.h>
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/nodes/NodeDisplacements.h"
#include "nuto/mechanics/loads/LoadNodeForces1D.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"

// constructor
NuTo::LoadNodeForces1D::LoadNodeForces1D(const NodeBase* rNode, double rDirection, double rValue) :
        LoadNode(rNode)
{
    // set direction
    if (fabs(rDirection) < 1e-14)
    {
        throw MechanicsException("[NuTo::LoadNodeForces1D::LoadNodeForces1D] direction vector has zero length.");
    }
    this->mDirection = rDirection/fabs(rDirection);

    // set value
    this->mValue = rValue;
}

// adds the load to global sub-vectors
void NuTo::LoadNodeForces1D::AddLoadToGlobalSubVectors(NuTo::FullMatrix<double>& rActiceDofsLoadVector, NuTo::FullMatrix<double>& rDependentDofsLoadVector)const
{
    assert(rActiceDofsLoadVector.GetNumColumns()==1);
    assert(rDependentDofsLoadVector.GetNumColumns()==1);
    try
    {
        const NodeDisplacements<1>& node = dynamic_cast<const NodeDisplacements<1> &>(*(this->mNode));
        int dof = node.GetDofDisplacement(0);
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
