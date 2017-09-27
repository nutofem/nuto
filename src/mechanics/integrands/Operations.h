#pragma once

#include "mechanics/nodes/DofVector.h"
#include "mechanics/nodes/DofMatrix.h"
#include "mechanics/integrands/Base.h"

namespace NuTo
{

class CellData;
class CellIpData;

struct VectorOperation
{
    virtual DofVector<double> operator()(Integrand::Base&, const CellData&, const CellIpData&) const = 0;
};

struct MatrixOperation
{
    virtual DofMatrix<double> operator()(Integrand::Base&, const CellData&, const CellIpData&) const = 0;
};

struct VoidOperation
{
    virtual void operator()(Integrand::Base&, const CellData&, const CellIpData&) const = 0;
};

} /* NuTo */
