#pragma once

#include "mechanics/constitutive/laws/LinearElastic.h"
#include "mechanics/dofs/DofType.h"
#include "mechanics/dofs/DofVector.h"
#include "mechanics/dofs/DofMatrix.h"
#include "mechanics/interpolation/TypeDefs.h"
#include "mechanics/cell/CellIpData.h"
#include "mechanics/cell/CellData.h"

namespace NuTo
{
namespace Integrands
{

template <int TDim>
class MomentumBalance
{
public:
    MomentumBalance(NuTo::DofType dofType, const Laws::MechanicsInterface<TDim>& law)
        : mDofType(dofType)
        , mLaw(law)
    {
    }

    NuTo::DofVector<double> Gradient(const NuTo::CellData& cellData, const NuTo::CellIpData& cellIpData, double deltaT)
    {
        NuTo::BMatrixStrain B = cellIpData.GetBMatrixStrain(mDofType);
        NuTo::NodeValues u = cellData.GetNodeValues(mDofType);
        NuTo::DofVector<double> gradient;

        NuTo::EngineeringStrainPDE<TDim> strain = B * u;
        gradient[mDofType] = B.transpose() * mLaw.Stress(strain, deltaT, cellData.GetCellId(), cellIpData.GetIPNum());
        return gradient;
    }

    NuTo::DofMatrix<double> Hessian0(const NuTo::CellData& cellData, const NuTo::CellIpData& cellIpData, double deltaT)
    {
        NuTo::BMatrixStrain B = cellIpData.GetBMatrixStrain(mDofType);
        NuTo::NodeValues u = cellData.GetNodeValues(mDofType);
        NuTo::DofMatrix<double> hessian0;

        NuTo::EngineeringStrainPDE<TDim> strain = B * u;
        hessian0(mDofType, mDofType) =
                B.transpose() * mLaw.Tangent(strain, deltaT, cellData.GetCellId(), cellIpData.GetIPNum()) * B;
        return hessian0;
    }

private:
    NuTo::DofType mDofType;
    const Laws::MechanicsInterface<TDim>& mLaw;
};
} /* Integrand */
} /* NuTo */
