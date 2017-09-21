#pragma once
#include "mechanics/cell/CellData.h"
#include "mechanics/cell/CellIpData.h"
#include "mechanics/nodes/DofVector.h"
#include "mechanics/nodes/DofMatrix.h"
#include "mechanics/cell/IPValue.h"

namespace NuTo
{
template <int TDim>
class Integrand
{
public:
    virtual Integrand<TDim>* Clone() const = 0;
    virtual ~Integrand() = default;
    virtual DofVector<double> Gradient(const CellData&, const CellIpData<TDim>&) = 0;
    virtual DofMatrix<double> Hessian0(const CellData&, const CellIpData<TDim>&) = 0;
    virtual std::vector<IPValue> IPValues(const CellData&, const CellIpData<TDim>&) = 0;
};

} /* NuTo */
