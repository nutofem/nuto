#pragma once
#include "mechanics/cell/CellData.h"
#include "mechanics/cell/CellIPData.h"
#include "mechanics/nodes/DofVector.h"
#include "mechanics/nodes/DofMatrix.h"
#include "mechanics/cell/IPValue.h"

namespace NuTo
{
template <int TDim>
class Integrand
{
public:
    virtual Integrand* Clone() const = 0;
    virtual ~Integrand() = default;
    virtual DofVector<double> Gradient(const CellData<TDim>&, const CellIPData<TDim>&) = 0;
    virtual DofMatrix<double> Hessian0(const CellData<TDim>&, const CellIPData<TDim>&) = 0;
    virtual std::vector<IPValue> IPValues(const CellData<TDim>&, const CellIPData<TDim>&) = 0;
};

} /* NuTo */