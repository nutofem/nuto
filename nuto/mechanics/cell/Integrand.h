#pragma once
#include "nuto/mechanics/cell/CellData.h"
#include "nuto/mechanics/cell/CellIpData.h"
#include "nuto/mechanics/dofs/DofVector.h"
#include "nuto/mechanics/dofs/DofMatrix.h"
#include "nuto/mechanics/cell/IPValue.h"

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
