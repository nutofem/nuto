#pragma once

#include <Eigen/Core>
#include "nuto/base/Exception.h"

namespace NuTo
{

//! @brief stores the support points
class SupportPoints
{
public:
    //! @brief get number of support points
    inline int GetNumSupportPoints() const
    {
        return mSPOrigInput.cols();
    }

    //! @brief get output dimension of support points (input)
    inline int GetDimOutput() const
    {
        return mSPOrigOutput.rows();
    }

    //! @brief returns the input of the support points in a matrix
    inline const Eigen::MatrixXd& GetOrigSupportPointsInput() const
    {
        return mSPOrigInput;
    }

    //! @brief returns the output of the support points in a matrix
    inline const Eigen::MatrixXd& GetOrigSupportPointsOutput() const
    {
        return mSPOrigOutput;
    }

    void SetSupportPoints(const Eigen::MatrixXd& SPOrigInput, const Eigen::MatrixXd& SPOrigOutput)
    {
        if (SPOrigInput.cols() != SPOrigOutput.cols())
            throw Exception(__PRETTY_FUNCTION__,
                            "Number of columns for input and output must be identical (=number of samples).");

        mSPOrigInput = SPOrigInput;
        mSPOrigOutput = SPOrigOutput;
    }

private:
    Eigen::MatrixXd mSPOrigInput; //!< original inputs, each sample after another
    Eigen::MatrixXd mSPOrigOutput; //!< original outputs, each sample after another
};
} // namespace nuto
