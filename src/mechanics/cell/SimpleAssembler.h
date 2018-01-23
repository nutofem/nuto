#pragma once

#include "base/Group.h"
#include "mechanics/cell/CellInterface.h"
#include "mechanics/dofs/DofContainer.h"
#include "mechanics/dofs/GlobalDofVector.h"
#include "mechanics/dofs/GlobalDofMatrixSparse.h"

namespace NuTo
{
class SimpleAssembler
{
public:
    SimpleAssembler(DofContainer<int> numIndependentDofs, DofContainer<int> numDependentDofs);

    GlobalDofVector BuildVector(const Group<CellInterface>& cells, std::vector<DofType> dofTypes,
                                CellInterface::VectorFunction f) const;

    GlobalDofMatrixSparse BuildMatrix(const Group<CellInterface>& cells, std::vector<DofType> dofTypes,
                                      CellInterface::MatrixFunction f) const;

private:
    GlobalDofVector ProperlyResizedGlobalVector(std::vector<DofType> dofTypes) const;
    GlobalDofMatrixSparse ProperlyResizedGlobalMatrix(std::vector<DofType> dofTypes) const;

    NuTo::DofContainer<int> mNumIndependentDofs;
    NuTo::DofContainer<int> mNumDependentDofs;
};
} /* NuTo */
