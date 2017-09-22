#pragma once

#include "mechanics/cell/IntegrandTimeDependent.h"
#include "mechanics/cell/MechanicsLaw.h"
#include "mechanics/constitutive/EngineeringStrainPDE.h"
#include "mechanics/constitutive/EngineeringStressPDE.h"

namespace NuTo
{


template <int TDim>
class IntegrandMomentumBalance : public NuTo::IntegrandTimeDependent<TDim>
{
public:
    IntegrandMomentumBalance(const NuTo::DofType& dofType, const NuTo::MechanicsLaw<TDim>& law)
        : mDofType(dofType)
        , mLaw(law.Clone())
    {
    }

    virtual NuTo::IntegrandTimeDependent<TDim>* Clone() const override
    {
        return new IntegrandMomentumBalance<TDim>(mDofType, mLaw);
    }

    NuTo::DofVector<double> Gradient(const NuTo::CellData<TDim>& cellData,
                                     const NuTo::CellIPData<TDim>& cellIpData) override
    {
        NuTo::BMatrixStrain B = cellIpData.GetBMatrixStrain(mDofType);
        NuTo::NodeValues u = cellData.ExtractNodeValues(mDofType);
        NuTo::DofVector<double> gradient;

        NuTo::EngineeringStrainPDE<TDim> strain = B * u;
        gradient[mDofType] = B.transpose() * mLaw.Stress(strain);
        return gradient;
    }
    NuTo::DofMatrix<double> Hessian0(const NuTo::CellData<TDim>& cellData,
                                     const NuTo::CellIPData<TDim>& cellIpData) const override
    {
        NuTo::BMatrixStrain B = cellIpData.GetBMatrixStrain(mDofType);
        NuTo::DofMatrix<double> hessian0;

        NuTo::EngineeringStrainPDE<TDim> dummy;
        hessian0(mDofType, mDofType) = B.transpose() * mLaw.Tangent(dummy) * B;
        return hessian0;
    }

//    std::vector<NuTo::IPValue> IPValues(double t, double delta_t,
//                                        const NuTo::CellData& cellData,
//                                        const NuTo::CellIPData& cellIpData) const override
//    {
//        NuTo::BMatrixStrain B = cellIpData.GetBMatrixStrain(mDofType);
//        NuTo::NodeValues u = cellData.ExtractNodeValues(mDofType);
//
//        std::vector<NuTo::IPValue> ipValues;
//        Eigen::VectorXd strain = B * u;
//        Eigen::VectorXd stress = mLaw.Stress(strain);
//
//        ipValues.push_back({"Stress", stress});
//        ipValues.push_back({"Strain", strain});
//
//        return ipValues;
//    }

private:
    const NuTo::DofType& mDofType;
    MechanicsLaw<TDim>& mLaw;
};

} /* NuTo */
