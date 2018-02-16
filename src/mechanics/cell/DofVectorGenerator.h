#pragma once
#include <vector>
#include "base/EntryGenerator.h"
#include "mechanics/dofs/DofMatrix.h"
#include "mechanics/dofs/DofVector.h"

namespace NuTo
{
using VectorEntry = std::pair<DofVector<double>, DofVector<int>>;

struct DofVectorGenerator : EntryGenerator<VectorEntry>
{
    virtual ~DofVectorGenerator() = default;
    virtual std::vector<DofType> Dofs() const = 0;
};
} /* NuTo */
