// $Id$

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/loads/LoadNodeGroupForces3D.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"

//! @brief constructor
NuTo::LoadNodeGroupForces3D::LoadNodeGroupForces3D(const Group<NodeBase>* rGroup, const NuTo::FullMatrix<double>& rDirection, double rValue) :
        LoadNodeGroup(rGroup)
{
    if (rDirection.GetNumColumns()!=1 || rDirection.GetNumRows()!=3)
        throw MechanicsException("[NuTo::LoadNodeGroupForces3D::LoadNodeGroupForces3D] Dimension of the direction matrix must be equal to the dimension of the structure.");

    memcpy(mDirection,rDirection.mEigenMatrix.data(),3*sizeof(double));
    //normalize the direction
    double norm = sqrt(mDirection[0]*mDirection[0]+mDirection[1]*mDirection[1]+mDirection[2]*mDirection[2]);
    if (norm < 1e-14)
    {
        throw MechanicsException("[NuTo::LoadNodeGroupForces3D::LoadNodeGroupForces3D] direction vector has zero length");
    }
    double invNorm = 1./norm;
    mDirection[0]*=invNorm;
    mDirection[1]*=invNorm;
    mDirection[2]*=invNorm;
    mValue = rValue;
}

// adds the load to global sub-vectors
void NuTo::LoadNodeGroupForces3D::AddLoadToGlobalSubVectors(NuTo::FullMatrix<double>& rActiceDofsLoadVector, NuTo::FullMatrix<double>& rDependentDofsLoadVector)const
{
    assert(rActiceDofsLoadVector.GetNumColumns()==1);
    assert(rDependentDofsLoadVector.GetNumColumns()==1);
    for (Group<NodeBase>::const_iterator itNode=this->mGroup->begin(); itNode!=this->mGroup->end(); itNode++)
    {
        try
        {
            for (int dofCount = 0; dofCount < 3; dofCount++)
            {
                int dof = (*itNode)->GetDofDisplacement(dofCount);
                assert(dof >= 0);
                if (dof < rActiceDofsLoadVector.GetNumRows())
                {
                    rActiceDofsLoadVector(dof,0) += this->mValue*mDirection[dofCount];
                }
                else
                {
                    dof -= rActiceDofsLoadVector.GetNumRows();
                    assert(dof < rDependentDofsLoadVector.GetNumRows());
                    rDependentDofsLoadVector(dof,0) += this->mValue*mDirection[dofCount];
                }
            }
        }
        catch (std::bad_cast & b)
        {
            throw MechanicsException("[NuTo::LoadNodeGroupForces3D::AddLoad] Node has no displacements or its dimension is not equivalent to the 1.");
        }
        catch (...)
        {
            throw MechanicsException("[NuTo::LoadNodeGroupForces3D::AddLoad] Error getting displacements of node (unspecified exception).");
        }
    }
}
