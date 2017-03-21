#pragma once

#include <eigen3/Eigen/Dense>
#include <cassert>

namespace NuTo
{
class NodeSimple
{
public:
    NodeSimple(Eigen::VectorXd rValues)
        : mValues(rValues)
    {
        mDofNumbers.setZero(mValues.rows());
    }

    const Eigen::VectorXd& GetValues() const
    {
        return mValues;
    }

    int GetDofNumber(int rComponent) const
    {
        assert(rComponent < mDofNumbers.rows());
        return mDofNumbers[rComponent];
    }

    void SetValue(int rComponent, double rValue)
    {
        assert(rComponent < mDofNumbers.rows());
        mValues[rComponent] = rValue;
    }

    void SetDofNumber(int rComponent, int rDofNumber)
    {
        assert(rComponent < mDofNumbers.rows());
        mDofNumbers[rComponent] = rDofNumber;
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
