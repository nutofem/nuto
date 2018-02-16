#pragma once

#include "mechanics/cell/DofVectorGenerator.h"
#include "mechanics/cell/DofMatrixGenerator.h"
#include "mechanics/dofs/DofInfo.h"
#include "mechanics/dofs/GlobalDofVector.h"
#include "mechanics/dofs/GlobalDofMatrixSparse.h"

namespace NuTo
{
class SimpleAssembler
{
public:
    SimpleAssembler() = default;
    SimpleAssembler(DofInfo dofInfo);

    GlobalDofVector BuildVector(DofVectorGenerator& entries) const;

    GlobalDofMatrixSparse BuildMatrix(DofMatrixGenerator& entries) const;

    void SetDofInfo(DofInfo dofInfo);

private:
    GlobalDofVector ProperlyResizedGlobalVector(std::vector<DofType> dofTypes) const;
    GlobalDofMatrixSparse ProperlyResizedGlobalMatrix(std::vector<DofType> dofTypes) const;

    DofInfo mDofInfo;

    void ThrowOnZeroDofNumbering(std::vector<DofType> dofTypes) const;
};
} /* NuTo */
