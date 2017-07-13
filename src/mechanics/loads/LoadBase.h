#pragma once

#include "mechanics/structures/StructureOutputBlockVector.h"

namespace NuTo
{

//! @brief Abstract class for all constraint equations
class LoadBase
{

public:
    //! @brief Constructor
    LoadBase()
    {
    }

    //! @brief Destructor
    virtual ~LoadBase()
    {
    }

    //! @brief Adds the load to global sub-vectors
    //! @param rActiceDofsLoadVector Global load vector which correspond to the active dofs
    //! @param rDependentDofsLoadVector Global load vector which correspond to the dependent dofs
    virtual void AddLoadToGlobalSubVectors(StructureOutputBlockVector& externalLoad) const = 0;
};
} // namespace NuTo
