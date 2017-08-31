#pragma once

#include "mechanics/cell/Integrand.h"
#include "mechanics/constitutive/laws/LinearElastic.h"

namespace NuTo
{


template <int TDim>
class IntegrandLinearElastic : public NuTo::Integrand<TDim>
{
public:
    IntegrandLinearElastic(const NuTo::DofType& rDofType, const Laws::LinearElastic<TDim>& rLaw)
        : mDofType(rDofType)
        , mLaw(rLaw)
    {
    }

    virtual NuTo::Integrand<TDim>* Clone() const override
    {
        return new IntegrandLinearElastic<TDim>(mDofType, mLaw);
    }

    NuTo::DofVector<double> Gradient(const NuTo::CellData& rCellData,
                                     const NuTo::CellIPData<TDim>& rCellIPData) override
    {
        NuTo::BMatrixStrain B = rCellIPData.GetBMatrixStrain(mDofType);
        NuTo::NodeValues u = rCellData.GetNodeValues(mDofType);
        NuTo::DofVector<double> gradient;

        NuTo::EngineeringStrainPDE<TDim> strain = B * u;
        gradient[mDofType] = B.transpose() * mLaw.Stress(strain);
        return gradient;
    }

    NuTo::DofMatrix<double> Hessian0(const NuTo::CellData& rCellData,
                                     const NuTo::CellIPData<TDim>& rCellIPData) override
    {
        NuTo::BMatrixStrain B = rCellIPData.GetBMatrixStrain(mDofType);
        NuTo::DofMatrix<double> hessian0;

        NuTo::EngineeringStrainPDE<TDim> dummy;
        hessian0(mDofType, mDofType) = B.transpose() * mLaw.Tangent(dummy) * B;
        return hessian0;
    }
    std::vector<NuTo::IPValue> IPValues(const NuTo::CellData& rCellData,
                                        const NuTo::CellIPData<TDim>& rCellIPData) override
    {
        NuTo::BMatrixStrain B = rCellIPData.GetBMatrixStrain(mDofType);
        NuTo::NodeValues u = rCellData.GetNodeValues(mDofType);

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
