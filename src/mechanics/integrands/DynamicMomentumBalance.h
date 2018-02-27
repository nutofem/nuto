#pragma once

#include "mechanics/integrands/MomentumBalance.h"

namespace NuTo
{
namespace Integrands
{

template <int TDim>
class DynamicMomentumBalance : public MomentumBalance<TDim>
{
public:
    DynamicMomentumBalance(DofType dofType, const Laws::MechanicsInterface<TDim>& law, double rho)
        : MomentumBalance<TDim>(dofType, law)
        , mRho(rho)
    {
    }

    DofMatrix<double> Hessian2(const CellIpData& cellIpData)
    {

        NMatrix N = cellIpData.N(MomentumBalance<TDim>::mDofType);
        DofMatrix<double> hessian2;

        hessian2(MomentumBalance<TDim>::mDofType, MomentumBalance<TDim>::mDofType) = N.transpose() * N * mRho;
        return hessian2;
    }

private:
    double mRho;
};
} /* Integrand */
} /* NuTo */
