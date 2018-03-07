#pragma once

#include <Eigen/Core>

namespace NuTo
{
namespace ShapeFunctions1D
{

//! @brief Lobatto nodes in one dimension
//!
//! Node Numbering:
//!
//! number               0   ,  1,  2, ... , order
//!                      |--------------------|
//! local coordinate    -1   ,  .......... , +1
//!
//! @param order Polynomial order (minumum: 1 linear, 2 quadratic, etc.)
//! @return local node coordinates ordered from left to right
Eigen::VectorXd NodeCoordinatesTrussLobatto(int order);


//! @brief Compute barycentric weights from node coordinates
//! @param nodes local node coordinates
Eigen::VectorXd BarycentricWeights(const Eigen::VectorXd& nodes);


//! @brief Value of 1D Lagrange type shape functions with given nodes
//!
//! Computed using barycentric formula
//! @param x local coordinate where shapes are evaluated
//! @param nodes local node coordinates
//! @return shapes evaluated at x. Fulfill interpolation condition fi(xj) = delta_ij
Eigen::VectorXd ShapeFunctionsTrussLagrange(const double x, const Eigen::VectorXd& nodes);


//! @brief Value of 1D Lagrange type shape function derivative with given nodes
//! @param x local coordinate where shapes are evaluated
//! @param nodes local node coordinates
//! @return shape derivatives evaluated at x.
Eigen::VectorXd DerivativeShapeFunctionsTrussLagrange(const double x, const Eigen::VectorXd& nodes);
}

namespace ShapeFunctions2D
{
//! @brief Lobatto nodes for quad element (tensor product)
//!
//! Node Numbering:
//!
//!                    y ^
//!                      |
//!
//!                     12---13------14------15
//!                      |                    |
//!                      |                    |
//!                      8    9      10      11
//!                      |                    |
//!                      |                    |
//!                      4     5      6       7
//!                      |                    |
//!                      |                    |
//!                      0-----1------2-------3     --> x
//!
//!
//! local coordinate    -1   ,  .......... , +1
//!
//! @param nodeId
//! @param nodes local node coords in interval [-1,1 ] used to build
//! a quad grid
//! @return coordinates of node with id nodeId
Eigen::MatrixXd NodeCoordinatesQuadLobatto(int nodeId, const Eigen::VectorXd& nodes);


//! @brief Value of Lagrange type shape functions with given 1D nodes (tensor product)
//! @param x local coordinate where shapes are evaluated
//! @param nodes local node coordinates
//! @return shapes evaluated at x. Fulfill interpolation condition fi(xj) = delta_ij
Eigen::VectorXd ShapeFunctionsQuadLagrange(const Eigen::Vector2d x, const Eigen::VectorXd& nodes);


Eigen::MatrixXd DerivativeShapeFunctionsQuadLagrange(const Eigen::Vector2d x, const Eigen::VectorXd& nodes);
}


namespace ShapeFunctions3D
{
//! @brief Lobatto nodes for brick element (tensor product)
//!
//! Node Numbering:
//!
//!                              z
//!                             /
//!
//!
//!
//!                    y ^  /                   /
//!                      | 32---33------34------35
//!                       /                   /
//!                     12---13------14------15
//!                      |                    | 27
//!                      |                    |
//!                      8    9      10      11
//!                      |                    | 23
//!                      |                    |
//!                      4     5      6       7  /
//!                      |                    | 19
//!                      |                    |/
//!                      0-----1------2-------3     --> x
//!
//!
//! local coordinate    -1   ,  .......... , +1
//!
//! @param nodeId
//! @param nodes local node coords in interval [-1,1 ] used to build
//! a brick grid
//! @return coordinates of node with id nodeId
Eigen::MatrixXd NodeCoordinatesBrickLobatto(int nodeId, const Eigen::VectorXd& nodes);


//! @brief Value of Lagrange type shape functions with given 1D nodes (tensor product)
//! @param x local coordinate where shapes are evaluated
//! @param nodes local node coordinates
//! @return shapes evaluated at x. Fulfill interpolation condition fi(xj) = delta_ij
Eigen::VectorXd ShapeFunctionsBrickLagrange(const Eigen::Vector3d x, const Eigen::VectorXd& nodes);


Eigen::MatrixXd DerivativeShapeFunctionsBrickLagrange(const Eigen::Vector3d x, const Eigen::VectorXd& nodes);
}

} /* namespace NuTo */
