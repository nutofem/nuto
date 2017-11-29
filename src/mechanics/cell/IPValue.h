#pragma once

#include <string>
#include <Eigen/Core>

namespace NuTo
{
struct IPValue
{
    std::string mName;
    Eigen::MatrixXd mValue;
};
} /* NuTo */
