#pragma once

#include "mechanics/cell/Integrand.h"

namespace NuTo
{


class LinearElasticLaw2D
{
public:
    LinearElasticLaw2D(double rE, double rNu)
        : mE(rE)
        , mNu(rNu)
    {
        double factor = mE / (1.0 - (mNu * mNu));
        double C11    = factor;
        double C12    = factor * mNu;
        double C33    = factor * 0.5 * (1.0 - mNu);

        mC = Eigen::Matrix3d::Zero();
        mC(0, 0) = C11;
        mC(1, 0) = C12;

        mC(0, 1) = C12;
        mC(1, 1) = C11;

        mC(2, 2) = C33;
    }

    Eigen::Vector3d Stress(Eigen::Vector3d rStrain) const
    {
        return mC * rStrain;
    }

    const Eigen::Matrix3d& C() const
    {
        return mC;
    }

    double mE;
    double mNu;
    Eigen::Matrix3d mC;
};

template <int TDim>
class IntegrandLinearElastic : public NuTo::Integrand<TDim>
{
public:
    IntegrandLinearElastic(const NuTo::DofType& rDofType, const LinearElasticLaw2D& rLaw)
        : mDofType(rDofType)
        , mLaw(rLaw)
    {
    }

    virtual NuTo::Integrand<TDim>* Clone() const override
    {
        return new IntegrandLinearElastic<TDim>(mDofType,mLaw);
    }

    NuTo::DofVector<double> Gradient(const NuTo::CellData& rCellData,
                                     const NuTo::CellIPData<TDim>& rCellIPData) override
    {
        NuTo::BMatrixStrain B = rCellIPData.GetBMatrixStrain(mDofType);
        NuTo::NodeValues u    = rCellData.GetNodeValues(mDofType);
        NuTo::DofVector<double> gradient;

        Eigen::VectorXd strain = B * u;
        gradient[mDofType]     = B.transpose() * mLaw.Stress(strain);
        return gradient;
    }

    std::vector<NuTo::IPValue> IPValues(const NuTo::CellData& rCellData,
                                        const NuTo::CellIPData<TDim>& rCellIPData) override
    {
        NuTo::BMatrixStrain B = rCellIPData.GetBMatrixStrain(mDofType);
        NuTo::NodeValues u    = rCellData.GetNodeValues(mDofType);

        std::vector<NuTo::IPValue> ipValues;
        Eigen::VectorXd strain = B * u;
        Eigen::VectorXd stress = mLaw.Stress(strain);

        ipValues.push_back({"Stress", stress});
        ipValues.push_back({"Strain", strain});

        return ipValues;
    }

private:
    const NuTo::DofType& mDofType;
    const LinearElasticLaw2D& mLaw;
};

} /* NuTo */
