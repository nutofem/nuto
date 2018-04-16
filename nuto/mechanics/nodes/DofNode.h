#pragma once

#include <Eigen/Core>
#include <vector>
#include <cassert>

namespace NuTo
{
//! @brief Store node values and its dof
//! @todo fix sized nodes?
class DofNode
{
public:
    //! @brief initizalizes the node values with `values` and initializes the dof numbers to zero
    //! @param values inititial node values
    //! @remark this magic number `-1` indicates an uninitialized state. I was not able to declare a static variable
    //! NOT_SET=-1. Maybe someone else can help.
    DofNode(Eigen::VectorXd values)
        : mValues({values})
        , mDofNumbers(Eigen::VectorXi::Constant(values.rows(), -1))
    {
    }

    //! @brief initializes a 1D node with `value` and a dof number 0
    //! @param value initial node value
    DofNode(double value)
        : mValues({Eigen::VectorXd::Constant(1, value)})
        , mDofNumbers(Eigen::VectorXi::Zero(1))
    {
    }

    DofNode(int dimension, int numInstances)
        : mValues(numInstances, Eigen::VectorXd::Zero(dimension))
        , mDofNumbers(Eigen::VectorXi::Constant(dimension, -1))
    {
    }

    //! Allocates `numInstances` full of zeros
    void AllocateInstances(int numInstances)
    {
        assert(numInstances > 0);
        mValues.resize(numInstances, Eigen::VectorXd::Zero(mValues[0].size()));
    }

    int GetNumInstances() const
    {
        return mValues.size();
    }

    const Eigen::VectorXd& GetValues(int instance = 0) const
    {
        assert(instance < static_cast<int>(mValues.size()));
        return mValues[instance];
    }

    int GetDofNumber(int component) const
    {
        assert(component < mDofNumbers.rows());
        return mDofNumbers[component];
    }

    void SetValues(Eigen::VectorXd values, int instance = 0)
    {
        assert(instance < static_cast<int>(mValues.size()));
        assert(values.size() == mValues[instance].size());
        mValues[instance] = values;
    }

    void SetValue(int component, double value, int instance = 0)
    {
        assert(instance < static_cast<int>(mValues.size()));
        assert(component < mDofNumbers.rows());
        mValues[instance][component] = value;
    }

    void SetDofNumber(int component, int dofNumber)
    {
        assert(component < mDofNumbers.rows());
        mDofNumbers[component] = dofNumber;
    }

    int GetNumValues() const
    {
        return mValues[0].rows();
    }

private:
    std::vector<Eigen::VectorXd> mValues;
    Eigen::VectorXi mDofNumbers;
};
} /* NuTo */
