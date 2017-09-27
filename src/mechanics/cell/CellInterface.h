#pragma once

#include "mechanics/integrands/Operations.h"
#include "mechanics/nodes/DofVector.h"
#include "mechanics/nodes/DofMatrix.h"

namespace NuTo
{
class CellInterface
{
public:
    virtual ~CellInterface() = default;

    virtual DofVector<double> operator()(const VectorOperation&) = 0;
    virtual DofMatrix<double> operator()(const MatrixOperation&) = 0;
    virtual double operator()(const ScalarOperation&) = 0;
    virtual void operator()(const VoidOperation&) = 0;
    virtual DofVector<int> DofNumbering() = 0;
};
} /* NuTo */
