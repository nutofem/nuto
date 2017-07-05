// $Id$
#pragma once

#include "mechanics/loads/LoadNode.h"
#include "mechanics/structures/StructureOutputBlockVector.h"

namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... class for all scalar source terms (e.g. charges) applied to a single node
class LoadNodeScalarSource : public LoadNode
{

public:
    //! @brief constructor
    //! @param rValue ... value of source
    LoadNodeScalarSource(const NodeBase* rNode, double rValue);

    //! @brief adds the load to global sub-vectors
    //! @param rActiceDofsLoadVector ... global load vector which correspond to the active dofs
    //! @param rDependentDofsLoadVector ... global load vector which correspond to the dependent dofs
    void AddLoadToGlobalSubVectors(StructureOutputBlockVector& externalLoad) const override;

protected:
    double mValue;     //!< prescribed value of the node

private:
    LoadNodeScalarSource(){}
};
}//namespace NuTo

