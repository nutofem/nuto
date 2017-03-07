#pragma once
#include <memory>
#include "mechanics/cell/CellData.h"
#include "mechanics/cell/CellIPData.h"
#include "mechanics/nodes/DofVector.h"
#include "mechanics/cell/IPValue.h"

namespace NuTo
{
template <int TDim>
class Integrand
{
public:
    virtual std::unique_ptr<Integrand> Clone() const = 0;
    virtual ~Integrand()        = default;
    virtual DofVector Gradient(const CellData&, const CellIPData<TDim>&) = 0;
    virtual std::vector<IPValue> IPValues(const CellData&, const CellIPData<TDim>&) = 0;
};
} /* NuTo */
