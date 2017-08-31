#pragma once

#include <eigen3/Eigen/Core>
#include <cassert>

namespace NuTo
{
class NodeSimple
{
public:
    NodeSimple(Eigen::VectorXd values)
        : mValues(values)
          , mDofNumbers(Eigen::VectorXd::Zero(values.rows()))
    {
    }

    const Eigen::VectorXd& GetValues() const
    {
        return mValues;
    }

    int GetDofNumber(int component) const
    {
        assert(component < mDofNumbers.rows());
        return mDofNumbers[component];
    }

    void SetValue(int component, double value)
    {
        assert(component < mDofNumbers.rows());
        mValues[component] = value;
    }

    void SetDofNumber(int component, int dofNumber)
    {
        assert(component < mDofNumbers.rows());
        mDofNumbers[component] = dofNumber;
    }

    int GetNumValues() const
    {
        return mValues.rows();
    }

private:
    Eigen::VectorXd mValues;
    Eigen::VectorXi mDofNumbers;
};
} /* NuTo */
