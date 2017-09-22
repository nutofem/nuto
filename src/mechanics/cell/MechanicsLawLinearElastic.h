#pragma once
#include "mechanics/cell/MechanicsLaw.h"
namespace NuTo
{
template <int TDim>
class MechanicsLawLinearElastic: public MechanicsLaw<TDim>
{
public:
    MechanicsLawLinearElastic(const std::array<double, 3>& rMaterialParameters)
        : mMaterialParameters(rMaterialParameters){};
    MechanicsLaw<TDim>& Clone() const override
    {
        return *(new MechanicsLawLinearElastic<TDim>(mMaterialParameters));
    }
    Eigen::VectorXd Stress(const Eigen::VectorXd& rStrain) override
    {
        return this->Tangent() * rStrain;
    }
    inline Eigen::MatrixXd Tangent(const Eigen::VectorXd& rStrain) const override
    {
        return Tangent();
    }
    Eigen::MatrixXd Tangent() const
    {
        Eigen::MatrixXd tangent(GetVoigtDim(TDim), GetVoigtDim(TDim));
        switch (TDim)
        {
            case 1:
                tangent << mMaterialParameters[0];
                break;
            case 2: // plane strain
            {
                double s = mMaterialParameters[0] / ((1. + mMaterialParameters[1]) * (1. - 2. * mMaterialParameters[1]));
                double C11 = s * (1. - mMaterialParameters[1]);
                double C12 = s * mMaterialParameters[1];
                double C44 = s * (1. - 2 * mMaterialParameters[1]) * 0.5;
                tangent << C11, C12, 0, C12, C11, 0, 0, 0, C44;
            }
                break;
            case 3:
            {
                double s = mMaterialParameters[0] / ((1. + mMaterialParameters[1]) * (1. - 2. * mMaterialParameters[1]));
                double C11 = s * (1. - mMaterialParameters[1]);
                double C12 = s * mMaterialParameters[1];
                double C44 = s * (1. - 2 * mMaterialParameters[1]) * 0.5;
                tangent << C11, C12, C12, 0, 0, 0, C12, C11, C12, 0, 0, 0, C12, C12, C11, 0, 0, 0, 0, 0, 0, C44, 0, 0, 0, 0,
                    0, 0, C44, 0, 0, 0, 0, 0, 0, C44;
            }
                break;
        }
        return tangent;
    }
    void UpdateHistoryData(const Eigen::VectorXd& rStrain) override
    {
    }
    double Density() const override
    {
        return mMaterialParameters[2];
    }
private:
    const std::array<double, 3>& mMaterialParameters; // E,nu,rho
};
}