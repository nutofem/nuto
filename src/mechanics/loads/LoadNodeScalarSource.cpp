// $Id$
#include <assert.h>
#include "mechanics/MechanicsException.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/loads/LoadNodeScalarSource.h"

#include "math/SparseMatrixCSRGeneral.h"

// constructor
NuTo::LoadNodeScalarSource::LoadNodeScalarSource(const NodeBase* rNode, double rValue) :
        LoadNode(rNode)
{
    this->mValue = rValue;
}

// adds the load to global sub-vectors
void NuTo::LoadNodeScalarSource::AddLoadToGlobalSubVectors(StructureOutputBlockVector& externalLoad)const
{
    assert(externalLoad.J[Node::eDof::ELECTRICPOTENTIAL].cols()==1);
    assert(externalLoad.K[Node::eDof::ELECTRICPOTENTIAL].cols()==1);
    try
    {
        int dof = mNode->GetDof(Node::eDof::ELECTRICPOTENTIAL, 0);
        assert(dof >= 0);
        if (dof < externalLoad.J[Node::eDof::ELECTRICPOTENTIAL].rows())
        {
            externalLoad.J[Node::eDof::ELECTRICPOTENTIAL](dof,0) += this->mValue;
        }
        else
        {
            dof -= externalLoad.J[Node::eDof::ELECTRICPOTENTIAL].rows();
            assert(dof < externalLoad.K[Node::eDof::ELECTRICPOTENTIAL].rows());
            externalLoad.K[Node::eDof::ELECTRICPOTENTIAL](dof,0) += this->mValue;
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

