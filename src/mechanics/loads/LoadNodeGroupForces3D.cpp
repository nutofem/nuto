// $Id$

#include "mechanics/MechanicsException.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/groups/Group.h"
#include "mechanics/loads/LoadNodeGroupForces3D.h"
#include "math/SparseMatrixCSRGeneral.h"

//! @brief constructor
NuTo::LoadNodeGroupForces3D::LoadNodeGroupForces3D(int rLoadCase, const Group<NodeBase>* rGroup, const Eigen::VectorXd& rDirection, double rValue) :
        LoadNodeGroup(rLoadCase,rGroup)
{
    if (rDirection.rows()!=3)
        throw MechanicsException("[NuTo::LoadNodeGroupForces3D::LoadNodeGroupForces3D] Dimension of the direction matrix must be equal to the dimension of the structure.");

    memcpy(mDirection,rDirection.data(),3*sizeof(double));
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
void NuTo::LoadNodeGroupForces3D::AddLoadToGlobalSubVectors(int rLoadCase, Eigen::VectorXd& rActiceDofsLoadVector, Eigen::VectorXd& rDependentDofsLoadVector)const
{
    if (rLoadCase!=mLoadCase)
    	return;
    assert(rActiceDofsLoadVector.cols()==1);
    assert(rDependentDofsLoadVector.cols()==1);
    for (Group<NodeBase>::const_iterator itNode=this->mGroup->begin(); itNode!=this->mGroup->end(); itNode++)
    {
        try
        {
            for (int dofCount = 0; dofCount < 3; dofCount++)
            {
                int dof = itNode->second->GetDof(Node::eDof::DISPLACEMENTS, dofCount);
                assert(dof >= 0);
                if (dof < rActiceDofsLoadVector.rows())
                {
                    rActiceDofsLoadVector(dof,0) += this->mValue*mDirection[dofCount];
                }
                else
                {
                    dof -= rActiceDofsLoadVector.rows();
                    assert(dof < rDependentDofsLoadVector.rows());
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
