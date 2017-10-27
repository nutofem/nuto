#pragma once

#include <vector>
#include <string>
#include <eigen3/Eigen/Core>

namespace NuTo
{

//! @brief single integration point (ip) value entry
struct IpValue
{
    Eigen::VectorXd data;
    std::string name;
};

//! @brief collection of ip values for the same ip
using IpValues = std::vector<IpValue>;

} /* NuTo */
