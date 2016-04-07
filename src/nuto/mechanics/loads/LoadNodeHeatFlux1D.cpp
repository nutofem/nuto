#include <assert.h>
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/loads/LoadNodeHeatFlux1D.h"
#include "nuto/math/FullMatrix.h"

// constructor
NuTo::LoadNodeHeatFlux1D::LoadNodeHeatFlux1D(int rLoadCase, const NodeBase* rNode,
        double rDirection, double rValue) : LoadNode(rLoadCase,rNode)
{
    // set direction
    if (fabs(rDirection) < 1e-14)
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Direction vector has zero length.");
    }
    this->mDirection = rDirection/fabs(rDirection);

    // set value
    this->mValue = rValue;
}

// adds the load to global sub-vectors
void NuTo::LoadNodeHeatFlux1D::AddLoadToGlobalSubVectors(int rLoadCase,
        NuTo::FullVector<double,Eigen::Dynamic>& rActiceDofsLoadVector,
        NuTo::FullVector<double,Eigen::Dynamic>& rDependentDofsLoadVector)const
{
    if (rLoadCase!=mLoadCase)
    	return;
    assert(rActiceDofsLoadVector.GetNumColumns()==1);
    assert(rDependentDofsLoadVector.GetNumColumns()==1);
    try
    {
        int dof = mNode->GetDofTemperature();
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
        throw MechanicsException(__PRETTY_FUNCTION__, "Node has no temperature or its dimension is not equivalent to the 1.");
    }
    catch (...)
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Error getting temperature of node (unspecified exception).");
    }
}
