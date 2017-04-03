// $Id$
#include <assert.h>
#include "mechanics/MechanicsException.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/loads/LoadNodeScalarSource.h"

#include "math/SparseMatrixCSRGeneral.h"

// constructor
NuTo::LoadNodeScalarSource::LoadNodeScalarSource(int rLoadCase, const NodeBase* rNode, double rValue) :
        LoadNode(rLoadCase,rNode)
{
    this->mValue = rValue;
}

// adds the load to global sub-vectors
void NuTo::LoadNodeScalarSource::AddLoadToGlobalSubVectors(int rLoadCase, Eigen::VectorXd& rActiveDofsLoadVector, Eigen::VectorXd& rDependentDofsLoadVector)const
{
    assert(rActiveDofsLoadVector.cols()==1);
    assert(rDependentDofsLoadVector.cols()==1);
    try
    {
        int dof = mNode->GetDof(Node::eDof::ELECTRICPOTENTIAL, 0);
        assert(dof >= 0);
        if (dof < rActiveDofsLoadVector.rows())
        {
            rActiveDofsLoadVector(dof,0) += this->mValue;
        }
        else
        {
            dof -= rActiveDofsLoadVector.rows();
            assert(dof < rDependentDofsLoadVector.rows());
            rDependentDofsLoadVector(dof,0) += this->mValue;
        }
    }
    catch (std::bad_cast & b)
    {
        throw MechanicsException("[NuTo::LoadNodeScalarSource::AddLoad] Node has no electric potential or its dimension is not equivalent to the 1.");
    }
    catch (...)
    {
        throw MechanicsException("[NuTo::LoadNodeScalarSource::AddLoad] Error getting potential of node (unspecified exception).");
    }
}

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::LoadNodeScalarSource)
BOOST_CLASS_TRACKING(NuTo::LoadNodeScalarSource, track_always)
#endif
