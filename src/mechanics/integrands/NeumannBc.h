#pragma once

#include "mechanics/integrands/TimeDependent.h"
#include "mechanics/dofs/DofType.h"
#include "mechanics/interpolation/TypeDefs.h"
#include "mechanics/cell/CellIpData.h"
#include "mechanics/cell/CellData.h"
#include <iostream>

namespace NuTo
{
namespace Integrands
{
namespace TimeDependent
{

//! @brief applies a neumann boundary condition
//! @tparam TDim global dimension, say, node coordinate dimension. The dimension of the interpolation has to be one
//! order lower.
template <int TDim>
class NeumannBc : public Interface
{
public:
    NeumannBc(const NuTo::DofType& dofType, Eigen::Matrix<double, TDim, 1> factor)
        : mDofType(dofType)
        , mFactor(factor)
    {
    }

    std::unique_ptr<Base> Clone() const override
    {
        return std::make_unique<NeumannBc<TDim>>(*this);
    }

    NuTo::DofVector<double> Gradient(const NuTo::CellData& cellData, const NuTo::CellIpData& cellIpData, double = 0,
                                     double = 0) override
    {
        NuTo::NMatrix N = cellIpData.GetNMatrix(mDofType);
        NuTo::DofVector<double> gradient;

        gradient[mDofType] = N.transpose() * mFactor;

        return gradient;
    }

    NuTo::DofMatrix<double> Hessian0(const NuTo::CellData&, const NuTo::CellIpData& cellIpData, double = 0,
                                     double = 0) override
    {
        const int size = cellIpData.GetNMatrix(mDofType).cols();
        DofMatrix<double> hessian0;
        hessian0(mDofType, mDofType) = Eigen::MatrixXd::Zero(size, size);
        return hessian0;
    }

private:
    const NuTo::DofType& mDofType;
    Eigen::Matrix<double, TDim, 1> mFactor;
};

} /* TimeDependent */
} /* Integrand */
} /* NuTo */
