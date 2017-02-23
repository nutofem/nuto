#pragma once

#include <eigen3/Eigen/Dense>
#include <assert.h>

namespace NuTo
{
class NodeSimple
{
public:
    NodeSimple(Eigen::VectorXd rValues) : mValues(rValues)
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

    int GetNumValues() const
    {
        return mValues.rows();
    }

private:
    Eigen::VectorXd mValues;
    Eigen::VectorXi mDofNumbers;
};
} /* NuTo */
