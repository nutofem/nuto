#pragma once

#include "nuto/mechanics/dofs/DofType.h"
#include "nuto/mechanics/dofs/DofMatrix.h"
#include "nuto/mechanics/dofs/DofVector.h"
#include "nuto/mechanics/interpolation/TypeDefs.h"
#include "nuto/mechanics/cell/CellIpData.h"
#include "nuto/mechanics/cell/CellData.h"
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
    NeumannBc(NuTo::DofType dofType, Eigen::Matrix<double, TDim, 1> factor)
        : mDofType(dofType)
        , mFactor(factor)
    {
    }

    NuTo::DofVector<double> ExternalLoad(const NuTo::CellIpData& cellIpData)
    {
        Eigen::MatrixXd N = cellIpData.N(mDofType);
        DofVector<double> gradient;

        gradient[mDofType] = N.transpose() * mFactor;

        return gradient;
    }

private:
    DofType mDofType;
    Eigen::Matrix<double, TDim, 1> mFactor;
};
} /* Integrand */
} /* NuTo */
