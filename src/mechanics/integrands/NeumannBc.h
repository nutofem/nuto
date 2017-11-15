#pragma once

#include "mechanics/dofs/DofType.h"
#include "mechanics/dofs/DofMatrix.h"
#include "mechanics/dofs/DofVector.h"
#include "mechanics/interpolation/TypeDefs.h"
#include "mechanics/cell/CellIpData.h"
#include "mechanics/cell/CellData.h"

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
    NeumannBc(NuTo::DofType dofType, Eigen::Matrix<double, TDim, 1> factor)
        : mDofType(dofType)
        , mFactor(factor)
    {
    }

    NuTo::DofVector<double> Gradient(const NuTo::CellData& cellData, const NuTo::CellIpData& cellIpData)
    {
        NuTo::NMatrix N = cellIpData.GetNMatrix(mDofType);
        NuTo::DofVector<double> gradient;

        gradient[mDofType] = N.transpose() * mFactor;
        return gradient;
    }

    NuTo::DofMatrix<double> Hessian0(const NuTo::CellData&, const NuTo::CellIpData& cellIpData)
    {
        const int size = cellIpData.GetNMatrix(mDofType).cols();
        DofMatrix<double> hessian0;
        hessian0(mDofType, mDofType) = Eigen::MatrixXd::Zero(size, size);
        return hessian0;
    }

private:
    NuTo::DofType mDofType;
    Eigen::Matrix<double, TDim, 1> mFactor;
};
} /* Integrand */
} /* NuTo */
