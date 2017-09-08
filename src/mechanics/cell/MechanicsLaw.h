#pragma once

#include <eigen3/Eigen/Core>

namespace NuTo
{

template <int TDim>
class MechanicsLaw
{
public:
    virtual MechanicsLaw<TDim>& Clone() const = 0;

    virtual Eigen::VectorXd Stress(const Eigen::VectorXd& rStrain)=0;

    virtual Eigen::MatrixXd Tangent (const Eigen::VectorXd& rStrain) const = 0;

    virtual void UpdateHistoryData(const Eigen::VectorXd& rStrain) =0;

    virtual double Density() const = 0;
};


constexpr int GetVoigtDim(int rDim)
{
    switch (rDim)
    {
        case 1:
            return 1;
        case 2:
            return 3;
        case 3:
            return 6;
        default:
            throw;
    }
}
}
