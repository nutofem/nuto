#pragma once

#include "mechanics/integrands/Base.h"
#include "mechanics/integrands/Operations.h"
#include "mechanics/nodes/DofVector.h"
#include "mechanics/nodes/DofMatrix.h"

namespace NuTo
{
namespace Integrands
{
namespace TimeDependent
{
struct Interface : Base
{
    virtual ~Interface() = default;
    virtual DofVector<double> Gradient(const CellData&, const CellIpData&, double t, double deltaT) = 0;
    virtual DofMatrix<double> Hessian0(const CellData&, const CellIpData&, double t, double deltaT) = 0;
};

struct Gradient : VectorOperation
{
    Gradient(double t = 0, double deltaT = 0)
        : mT(t)
        , mDeltaT(deltaT)
    {
    }
    DofVector<double> operator()(Base& base, const CellData& cellData, const CellIpData& cellIpData) const override
    {
        return base.As<Interface>().Gradient(cellData, cellIpData, mT, mDeltaT);
    }

private:
    double mT;
    double mDeltaT;
};

struct Hessian0 : MatrixOperation
{
    Hessian0(double t = 0, double deltaT = 0)
        : mT(t)
        , mDeltaT(deltaT)
    {
    }
    DofMatrix<double> operator()(Base& base, const CellData& cellData, const CellIpData& cellIpData) const override
    {
        return base.As<Interface>().Hessian0(cellData, cellIpData, mT, mDeltaT);
    }
private:
    double mT;
    double mDeltaT;
};

} /*TimeDependent */
} /* Integrand */
} /* NuTo */
