#pragma once

#include <eigen3/Eigen/Core>

namespace NuTo
{

//! @brief Collection of helper functions for Eigen matrices
class EigenCompanion
{
public:
    //! @brief Append the rows of one matrix to the other
    //! @param top Matrix that is getting extended
    //! @param bottom Matrix that is added to the bottom of the other one
    static void AppendRows(Eigen::MatrixXd& top, const Eigen::MatrixXd& bottom);
};

}
