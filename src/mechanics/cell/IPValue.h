#pragma once

#include <string>
#include <eigen3/Eigen/Core>

namespace NuTo
{
struct IPValue
{
    std::string mName;
    Eigen::MatrixXd mValue;
};
} /* NuTo */
