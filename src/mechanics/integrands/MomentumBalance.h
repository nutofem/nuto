#pragma once

#include "mechanics/integrands/TimeDependent.h"
#include "mechanics/constitutive/laws/LinearElastic.h"
#include "mechanics/nodes/DofType.h"
#include "mechanics/interpolation/TypeDefs.h"
#include "mechanics/cell/CellIpData.h"
#include "mechanics/cell/CellData.h"

namespace NuTo
{
namespace Integrand
{
namespace TimeDependent
{

template <int TDim>
class MomentumBalance : public Interface
{
public:
    MomentumBalance(const NuTo::DofType& dofType, const Laws::MechanicsInterface<TDim>& law)
        : mDofType(dofType)
        , mLaw(law)
    {
    }

    std::unique_ptr<Base> Clone() const override
    {
        return std::make_unique<MomentumBalance<TDim>>(*this);
    }

    NuTo::DofVector<double> Gradient(const NuTo::CellData& cellData, const NuTo::CellIpData& cellIpData, double,
                                     double deltaT) override
    {
        NuTo::BMatrixStrain B = cellIpData.GetBMatrixStrain(mDofType);
        NuTo::NodeValues u = cellData.GetNodeValues(mDofType);
        NuTo::DofVector<double> gradient;

        NuTo::EngineeringStrainPDE<TDim> strain = B * u;
        gradient[mDofType] = B.transpose() * mLaw.Stress(strain, deltaT);
        return gradient;
    }

    NuTo::DofMatrix<double> Hessian0(const NuTo::CellData& cellData, const NuTo::CellIpData& cellIpData, double,
                                     double deltaT) override
    {
        NuTo::BMatrixStrain B = cellIpData.GetBMatrixStrain(mDofType);
        NuTo::NodeValues u = cellData.GetNodeValues(mDofType);
        NuTo::DofMatrix<double> hessian0;

        NuTo::EngineeringStrainPDE<TDim> strain = B * u;
        hessian0(mDofType, mDofType) = B.transpose() * mLaw.Tangent(strain, deltaT) * B;
        return hessian0;
    }

private:
    const NuTo::DofType& mDofType;
    const Laws::MechanicsInterface<TDim>& mLaw;
};

} /* TimeDependent */
} /* Integrand */
} /* NuTo */
