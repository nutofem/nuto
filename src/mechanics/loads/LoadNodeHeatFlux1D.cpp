#include <assert.h>
#include "mechanics/MechanicsException.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/loads/LoadNodeHeatFlux1D.h"
#include "math/FullMatrix.h"

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
        Eigen::VectorXd& rActiceDofsLoadVector,
        Eigen::VectorXd& rDependentDofsLoadVector)const
{
    if (rLoadCase!=mLoadCase)
    	return;
    assert(rActiceDofsLoadVector.cols()==1);
    assert(rDependentDofsLoadVector.cols()==1);
    try
    {
        int dof = mNode->GetDof(Node::eDof::TEMPERATURE);
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
        throw MechanicsException(__PRETTY_FUNCTION__, "Node has no temperature or its dimension is not equivalent to the 1.");
    }
    catch (...)
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Error getting temperature of node (unspecified exception).");
    }
}
