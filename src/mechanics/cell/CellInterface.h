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

    virtual DofVector<double> Integrate(const VectorOperation&) = 0;
    virtual DofMatrix<double> Integrate(const MatrixOperation&) = 0;
    virtual double Integrate(const ScalarOperation&) = 0;
    virtual void Apply(const VoidOperation&) = 0;
    virtual DofVector<int> DofNumbering() = 0;
};
} /* NuTo */
