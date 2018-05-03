#pragma once

#include "nuto/mechanics/integrands/MomentumBalance.h"

namespace NuTo
{
namespace Integrands
{

template <int TDim>
class DynamicMomentumBalance : public MomentumBalance<TDim>
{
using MomentumBalance<TDim>::mDofType;

public:
    DynamicMomentumBalance(DofType dofType, const Laws::MechanicsInterface<TDim>& law, double rho)
        : MomentumBalance<TDim>(dofType, law)
        , mRho(rho)
    {
    }

    DofMatrix<double> Hessian2(const CellIpData& cellIpData)
    {

        Eigen::MatrixXd N = cellIpData.N(mDofType);
        DofMatrix<double> hessian2;

        hessian2(mDofType, mDofType) = N.transpose() * N * mRho;
        return hessian2;
    }

private:
    double mRho;
};
} /* Integrand */
} /* NuTo */
