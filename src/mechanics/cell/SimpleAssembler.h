#pragma once

#include "base/Group.h"
#include "mechanics/cell/CellInterface.h"
#include "mechanics/dofs/DofNumbering.h"
#include "mechanics/dofs/GlobalDofVector.h"
#include "mechanics/dofs/GlobalDofMatrixSparse.h"

namespace NuTo
{
class SimpleAssembler
{
public:
    SimpleAssembler() = default;
    SimpleAssembler(DofInfo dofInfo);

    GlobalDofVector BuildVector(const Group<CellInterface>& cells, std::vector<DofType> dofTypes,
                                CellInterface::VectorFunction f) const;

    GlobalDofMatrixSparse BuildMatrix(const Group<CellInterface>& cells, std::vector<DofType> dofTypes,
                                      CellInterface::MatrixFunction f) const;

    void SetDofInfo(DofInfo dofInfo);

private:
    GlobalDofVector ProperlyResizedGlobalVector(std::vector<DofType> dofTypes) const;
    GlobalDofMatrixSparse ProperlyResizedGlobalMatrix(std::vector<DofType> dofTypes) const;

    DofInfo mDofInfo;

    void ThrowOnZeroDofNumbering(std::vector<DofType> dofTypes) const;
};
} /* NuTo */
