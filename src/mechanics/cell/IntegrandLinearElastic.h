#pragma once

#include "mechanics/cell/Integrand.h"
#include "mechanics/constitutive/laws/LinearElastic.h"

namespace NuTo
{


template <int TDim>
class IntegrandLinearElastic : public NuTo::Integrand<TDim>
{
public:
    IntegrandLinearElastic(const NuTo::DofType& dofType, const Laws::LinearElastic<TDim>& law)
        : mDofType(dofType)
        , mLaw(law)
    {
    }

    virtual NuTo::Integrand<TDim>* Clone() const override
    {
        return new IntegrandLinearElastic<TDim>(mDofType, mLaw);
    }

    NuTo::DofVector<double> Gradient(const NuTo::CellData& cellData,
                                     const NuTo::CellIpData<TDim>& cellIpData) override
    {
        NuTo::BMatrixStrain B = cellIpData.GetBMatrixStrain(mDofType);
        NuTo::NodeValues u = cellData.GetNodeValues(mDofType);
        NuTo::DofVector<double> gradient;

        NuTo::EngineeringStrainPDE<TDim> strain = B * u;
        gradient[mDofType] = B.transpose() * mLaw.Stress(strain);
        return gradient;
    }

    NuTo::DofMatrix<double> Hessian0(const NuTo::CellData& cellData,
                                     const NuTo::CellIpData<TDim>& cellIpData) override
    {
        NuTo::BMatrixStrain B = cellIpData.GetBMatrixStrain(mDofType);
        NuTo::DofMatrix<double> hessian0;

        NuTo::EngineeringStrainPDE<TDim> dummy;
        hessian0(mDofType, mDofType) = B.transpose() * mLaw.Tangent(dummy) * B;
        return hessian0;
    }
    std::vector<NuTo::IPValue> IPValues(const NuTo::CellData& cellData,
                                        const NuTo::CellIpData<TDim>& cellIpData) override
    {
        NuTo::BMatrixStrain B = cellIpData.GetBMatrixStrain(mDofType);
        NuTo::NodeValues u = cellData.GetNodeValues(mDofType);

        std::vector<NuTo::IPValue> ipValues;
        Eigen::VectorXd strain = B * u;
        Eigen::VectorXd stress = mLaw.Stress(strain);

        ipValues.push_back({"Stress", stress});
        ipValues.push_back({"Strain", strain});

        return ipValues;
    }

private:
    const NuTo::DofType& mDofType;
    const Laws::LinearElastic<TDim>& mLaw;
};

} /* NuTo */
