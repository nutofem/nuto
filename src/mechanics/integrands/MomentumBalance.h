#pragma once

#include "mechanics/constitutive/MechanicsInterface.h"
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
    MomentumBalance(DofType dofType, const Laws::MechanicsInterface<TDim>& law)
        : mDofType(dofType)
        , mLaw(law)
    {
    }

    DofVector<double> Gradient(const CellData& cellData, const CellIpData& cellIpData, double deltaT)
    {
        BMatrixStrain B = cellIpData.GetBMatrixStrain(mDofType);
        NodeValues u = cellData.GetNodeValues(mDofType);
        DofVector<double> gradient;

        NuTo::EngineeringStrain<TDim> strain = B * u;
        gradient[mDofType] = B.transpose() * mLaw.Stress(strain, deltaT, cellData.GetCellId(), cellIpData.GetIpId());
        return gradient;
    }

    DofMatrix<double> Hessian0(const CellData& cellData, const CellIpData& cellIpData, double deltaT)
    {
        BMatrixStrain B = cellIpData.GetBMatrixStrain(mDofType);
        NodeValues u = cellData.GetNodeValues(mDofType);
        DofMatrix<double> hessian0;

        NuTo::EngineeringStrain<TDim> strain = B * u;
        hessian0(mDofType, mDofType) =
                B.transpose() * mLaw.Tangent(strain, deltaT, cellData.GetCellId(), cellIpData.GetIpId()) * B;
        return hessian0;
    }

private:
    DofType mDofType;
    const Laws::MechanicsInterface<TDim>& mLaw;
};
} /* Integrand */
} /* NuTo */
