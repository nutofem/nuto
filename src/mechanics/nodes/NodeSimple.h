#pragma once

#include <eigen3/Eigen/Core>
#include <cassert>

namespace NuTo
{
//! @brief Store node values and its dof
//! @todo time derivatives? fix sized nodes?
class NodeSimple
{
public:

    //! @brief initizalizes the node values with \p values and initializes the dof numbers to zero
    //! @param values inititial node values
    NodeSimple(Eigen::VectorXd values)
        : mValues(values)
          , mDofNumbers(Eigen::VectorXi::Zero(values.rows()))
    {
    }
   
    //! @brief initializes a 1D node with \p value and a dof number 0
    //! @param value initial node value
    NodeSimple(double value)
        : mValues(Eigen::VectorXd::Constant(1, value))
          , mDofNumbers(Eigen::VectorXi::Zero(1))
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
