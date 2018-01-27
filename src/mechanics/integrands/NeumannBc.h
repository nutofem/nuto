#pragma once

#include "mechanics/dofs/DofType.h"
#include "mechanics/dofs/DofMatrix.h"
#include "mechanics/dofs/DofVector.h"
#include "mechanics/interpolation/TypeDefs.h"
#include "mechanics/cell/CellIpData.h"
#include "mechanics/cell/CellData.h"
#include <iostream>

namespace NuTo
{
namespace Integrands
{

//! @brief applies a neumann boundary condition
//! @tparam TDim global dimension, say, node coordinate dimension. The dimension of the interpolation has to be one
//! order lower.
template <int TDim>
class NeumannBc
{
public:
    using LoadFunction = std::function<Eigen::Matrix<double, TDim, 1>(Eigen::VectorXd)>;
    static LoadFunction Constant(Eigen::Matrix<double, TDim, 1> factor)
    {
        return [factor](Eigen::VectorXd) { return factor; };
    }

    NeumannBc(NuTo::DofType dofType, LoadFunction f)
        : mDofType(dofType)
        , mLoadFunction(f)
    {
    }

    NuTo::DofVector<double> ExternalLoad(const NuTo::CellData&, const NuTo::CellIpData& cellIpData)
    {
        NuTo::NMatrix N = cellIpData.GetNMatrix(mDofType);
        NuTo::DofVector<double> gradient;

        gradient[mDofType] = N.transpose() * mLoadFunction(cellIpData.GlobalCoordinates());

        return gradient;
    }

private:
    NuTo::DofType mDofType;
    LoadFunction mLoadFunction;
};
} /* Integrand */
} /* NuTo */
