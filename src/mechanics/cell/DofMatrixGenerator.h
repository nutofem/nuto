#pragma once
#include <vector>
#include "base/EntryGenerator.h"
#include "mechanics/dofs/DofMatrix.h"
#include "mechanics/dofs/DofVector.h"

namespace NuTo
{
using MatrixEntry = std::pair<DofMatrix<double>, DofVector<int>>;

struct DofMatrixGenerator : EntryGenerator<MatrixEntry>
{
    virtual ~DofMatrixGenerator() = default;
    virtual std::vector<DofType> Dofs() const = 0;
};
} /* NuTo */
