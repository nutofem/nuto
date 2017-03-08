#pragma once

#include <vector>
#include "mechanics/nodes/DofVector.h"
#include "mechanics/cell/IPValue.h"

namespace NuTo
{
class CellInterface
{
public:
    CellInterface()                     = default;
    virtual ~CellInterface()            = default;
    CellInterface(const CellInterface&) = default;
    CellInterface(CellInterface&&)      = default;
    CellInterface& operator=(const CellInterface&) = default;
    CellInterface& operator=(CellInterface&&) = default;

    virtual DofVector<double> Gradient() = 0;
    virtual DofVector<int> DofNumbering() = 0;
    virtual std::vector<std::vector<IPValue>> IPValues() = 0;
};
} /* NuTo */
