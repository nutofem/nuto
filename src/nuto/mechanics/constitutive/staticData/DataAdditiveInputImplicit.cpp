#include "DataAdditiveInputImplicit.h"

NuTo::Constitutive::StaticData::DataAdditiveInputImplicit::DataAdditiveInputImplicit()
    : mLocalInputs(),
      mStress(),
      mTime(0)
{
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

double &NuTo::Constitutive::StaticData::DataAdditiveInputImplicit::GetTime()
{
    return mTime;
}

double NuTo::Constitutive::StaticData::DataAdditiveInputImplicit::GetTime() const
{
    return mTime;
}
