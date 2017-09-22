#pragma once

#include <vector>
#include "mechanics/dofs/DofVector.h"
#include "mechanics/dofs/DofMatrix.h"
#include "mechanics/cell/IPValue.h"

namespace NuTo
{
class CellInterface
{
public:
    CellInterface() = default;
    virtual ~CellInterface() = default;
    CellInterface(const CellInterface&) = default;
    CellInterface(CellInterface&&) = default;
    CellInterface& operator=(const CellInterface&) = default;
    CellInterface& operator=(CellInterface&&) = default;

    virtual DofVector<double> Gradient() = 0;
    virtual DofMatrix<double> Hessian0() = 0;
    virtual DofVector<int> DofNumbering() = 0;
    virtual std::vector<std::vector<IPValue>> IPValues() = 0;
};
} /* NuTo */
