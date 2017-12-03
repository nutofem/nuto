#pragma once

#include <Eigen/Core>

namespace NuTo
{
namespace Constitutive
{
namespace StaticData
{

//! @brief Storing moisture transport static data.
class DataAdditiveInputImplicit
{
public:
    DataAdditiveInputImplicit();
    ~DataAdditiveInputImplicit();


    Eigen::VectorXd& GetStress();
    const Eigen::VectorXd& GetStress() const;

    Eigen::VectorXd& GetLocalInputs();
    const Eigen::VectorXd& GetLocalInputs() const;

    Eigen::VectorXd& GetLocalInputRates();
    const Eigen::VectorXd& GetLocalInputRates() const;

    double& GetTime();
    double GetTime() const;

private:
    Eigen::VectorXd mLocalInputs;
    Eigen::VectorXd mLocalInputRates;
    Eigen::VectorXd mStress;
    double mTime;
};

} // namespace StaticData
} // namespace Constitutive
} // namespace NuTo
