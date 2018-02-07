#pragma once

#include <vector>
#include <string>
#include <Eigen/Core>

namespace NuTo
{

//! @brief Collection of helper functions for Eigen matrices
namespace EigenCompanion
{
//! @brief converts data to a 3D vector, fills with zeros if needed
//! @param data vector of arbitrary size
//! @return 3D vector
inline Eigen::Vector3d To3D(const Eigen::VectorXd& data)
{
    const int dimension = data.rows();
    Eigen::Vector3d vector3d = Eigen::Vector3d::Zero();
    vector3d.block(0, 0, dimension, 1) = data;
    return vector3d;
}

//! @brief transforms an initializer_list to an Eigen::VectorXd
inline Eigen::VectorXd ToEigen(std::initializer_list<double> l)
{
    Eigen::VectorXd v(l.size());

    int position = 0;
    for (auto value : l)
        v[position++] = value;

    return v;
}

//! @brief transforms a single double number to an Eigen::VectorXd
inline Eigen::VectorXd ToEigen(double d)
{
    return Eigen::VectorXd::Constant(1, d);
}
} /* EigenCompanion */
} /* NuTo */
