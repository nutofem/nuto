#include "DataAdditiveInputImplicit.h"

NuTo::Constitutive::StaticData::DataAdditiveInputImplicit::DataAdditiveInputImplicit()
    : mLocalInputs(),
      mStress(),
      mTime(0)
{
    // https://forum.kde.org/viewtopic.php?f=74&t=118922  --- possible solutions to ambiguous call problems when using constructor like this: VectorXd(0)
    assert(mLocalInputs.size() == 0);
    assert(mStress.size() == 0);
}

NuTo::Constitutive::StaticData::DataAdditiveInputImplicit::~DataAdditiveInputImplicit()
{}

Eigen::VectorXd &NuTo::Constitutive::StaticData::DataAdditiveInputImplicit::GetStress()
{
    return mStress;
}

const Eigen::VectorXd &NuTo::Constitutive::StaticData::DataAdditiveInputImplicit::GetStress() const
{
    return mStress;
}

Eigen::VectorXd &NuTo::Constitutive::StaticData::DataAdditiveInputImplicit::GetLocalInputs()
{
    return mLocalInputs;
}

const Eigen::VectorXd &NuTo::Constitutive::StaticData::DataAdditiveInputImplicit::GetLocalInputs() const
{
    return mLocalInputs;
}

Eigen::VectorXd &NuTo::Constitutive::StaticData::DataAdditiveInputImplicit::GetLocalInputRates()
{
    return mLocalInputRates;
}

const Eigen::VectorXd &NuTo::Constitutive::StaticData::DataAdditiveInputImplicit::GetLocalInputRates() const
{
    return mLocalInputRates;
}

double &NuTo::Constitutive::StaticData::DataAdditiveInputImplicit::GetTime()
{
    return mTime;
}

double NuTo::Constitutive::StaticData::DataAdditiveInputImplicit::GetTime() const
{
    return mTime;
}
