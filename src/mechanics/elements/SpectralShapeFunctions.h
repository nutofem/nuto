#pragma once

#include "math/Legendre.h"
#include "base/Exception.h"
#include <Eigen/Core>

namespace NuTo
{
namespace ShapeFunctions1D
{

//! @brief Lobatto nodes in one dimension
//! @param order Polynomial order (minumum: 1 linear, 2 quadratic, etc.)
//! @return local node coordinates ordered from left to right
std::vector<double> NodeCoordinatesTrussLobatto(int order)
{
    if (order < 1)
        throw NuTo::Exception(__PRETTY_FUNCTION__, "Order too low. Must be 1 or higher");

    std::vector<double> points = NuTo::Math::Polynomial::LegendreDerivRoots(order);
    points.insert(points.begin(), -1.);
    points.push_back(1.);
    return points;
}

//! @brief Compute barycentric weights from node coordinates
//! @param nodes local node coordinates
std::vector<double> BarycentricWeights(const std::vector<double>& nodes)
{
    std::vector<double> w(nodes.size(), 1.);
    for (size_t j = 0; j < nodes.size(); j++)
        for (size_t i = 0; i < nodes.size(); i++)
        {
            if (i != j)
            {
                w[j] /= (nodes[j] - nodes[i]);
            }
        }
    return w;
}

//! @brief Value of 1D Lagrange type shape functions with given nodes
//!
//! Computed using barycentric formula
//! @param x local coordinate where shapes are evaluated
//! @param nodes local node coordinates
//! @return shapes evaluated at x. Fulfill interpolation condition fi(xj) = delta_ij
std::vector<double> ShapeFunctionsTrussLagrange(const double x, const std::vector<double>& nodes)
{
    std::vector<double> result(nodes.size(), 0.);
    std::vector<double> w = BarycentricWeights(nodes);
    // Check if x is a node (or near)
    bool xMatchesNode = false;
    for (size_t j = 0; j < nodes.size(); j++)
    {
        if (std::abs(x - nodes[j]) < 1.e-15)
        {
            result[j] = 1.;
            xMatchesNode = true;
        }
    }
    if (xMatchesNode)
    {
        return result;
    }

    double sum = 0;

    for (size_t j = 0; j < nodes.size(); j++)
    {
        double tmp = w[j]/(x-nodes[j]);
        result[j] = tmp;
        sum += tmp;
    }

    for (size_t j = 0; j < nodes.size(); j++)
    {
        result[j] /= sum;
    }
    return result;
}

//! @brief Value of 1D Lagrange type shape function derivative with given nodes
//! @param x local coordinate where shapes are evaluated
//! @param nodes local node coordinates
//! @return shape derivatives evaluated at x.
std::vector<double> DerivativeShapeFunctionsTrussLagrange(const double x, const std::vector<double>& nodes)
{
    std::vector<double> result(nodes.size(), 0.);
    std::vector<double> w = BarycentricWeights(nodes);
    for (size_t j = 0; j < nodes.size(); j++)
    {
        for (size_t k = 0; k < nodes.size(); k++)
        {
            if (k != j)
            {
                double tmp = 1.;
                for (size_t i = 0; i < nodes.size(); i++)
                {
                    if ((i != j) && (i != k))
                    {
                        tmp *= (x - nodes[i]);
                    }
                }
                result[j] += tmp;
            }
        }
        result[j] *= w[j];
    }
    return result;
}
}

namespace ShapeFunctions2D
{
//! @brief Lobatto nodes for quad element (tensor product)
//! @param order Polynomial order in 1D
//! @return local node coordinates
std::vector<Eigen::Vector2d> NodeCoordinatesQuadLobatto(int order)
{
    std::vector<double> nodes = ShapeFunctions1D::NodeCoordinatesTrussLobatto(order);
    std::vector<Eigen::Vector2d> result;
    for (size_t i=0; i<nodes.size(); i++)
    {
        for (size_t j=0; j<nodes.size(); j++)
        {
            result.push_back(Eigen::Vector2d({nodes[i],nodes[j]}));
        }
    }
    return result;
}

//! @brief Value of Lagrange type shape functions with given 1D nodes (tensor product)
//! @param x local coordinate where shapes are evaluated
//! @param nodes local node coordinates
//! @return shapes evaluated at x. Fulfill interpolation condition fi(xj) = delta_ij
std::vector<double> ShapeFunctionsQuadLagrange(const Eigen::Vector2d x, const std::vector<double>& nodes)
{
    std::vector<double> Nx = ShapeFunctions1D::ShapeFunctionsTrussLagrange(x(0), nodes);
    std::vector<double> Ny = ShapeFunctions1D::ShapeFunctionsTrussLagrange(x(1), nodes);

    std::vector<double> result;
    for (size_t j=0; j<nodes.size(); j++)
    {
        for (size_t i=0; i<nodes.size(); i++)
        {
            result.push_back(Nx[i] * Ny[j]);
        }
    }
    return result;
}

std::vector<Eigen::Vector2d> DerivativeShapeFunctionsQuadLagrange(const Eigen::Vector2d x, const std::vector<double>& nodes)
{
    std::vector<double> DNx = ShapeFunctions1D::DerivativeShapeFunctionsTrussLagrange(x(0), nodes);
    std::vector<double> DNy = ShapeFunctions1D::DerivativeShapeFunctionsTrussLagrange(x(1), nodes);

    std::vector<double> Nx = ShapeFunctions1D::ShapeFunctionsTrussLagrange(x(0), nodes);
    std::vector<double> Ny = ShapeFunctions1D::ShapeFunctionsTrussLagrange(x(1), nodes);

    std::vector<Eigen::Vector2d> result;
    for (size_t j=0; j<nodes.size(); j++)
    {
        for (size_t i=0; i<nodes.size(); i++)
        {
            result.push_back( Eigen::Vector2d({DNx[i]*Ny[j] , Nx[i]*DNy[j]})   );
        }
    }
    return result;
}

}



namespace ShapeFunctions3D
{
//! @brief Lobatto nodes for brick element (tensor product)
//! @param order Polynomial order in 1D
//! @return local node coordinates
std::vector<Eigen::Vector3d> NodeCoordinatesBrickLobatto(int order)
{
    std::vector<double> nodes = ShapeFunctions1D::NodeCoordinatesTrussLobatto(order);
    std::vector<Eigen::Vector3d> result;
    for (size_t i=0; i<nodes.size(); i++)
    {
        for (size_t j=0; j<nodes.size(); j++)
        {
            for (size_t k=0; k<nodes.size(); k++)
            {
                result.push_back(Eigen::Vector3d({nodes[i],nodes[j],nodes[k]}));
            }
        }
    }
    return result;
}

//! @brief Value of Lagrange type shape functions with given 1D nodes (tensor product)
//! @param x local coordinate where shapes are evaluated
//! @param nodes local node coordinates
//! @return shapes evaluated at x. Fulfill interpolation condition fi(xj) = delta_ij
std::vector<double> ShapeFunctionsBrickLagrange(const Eigen::Vector3d x, const std::vector<double>& nodes)
{
    std::vector<double> Nx = ShapeFunctions1D::ShapeFunctionsTrussLagrange(x(0), nodes);
    std::vector<double> Ny = ShapeFunctions1D::ShapeFunctionsTrussLagrange(x(1), nodes);
    std::vector<double> Nz = ShapeFunctions1D::ShapeFunctionsTrussLagrange(x(2), nodes);

    std::vector<double> result;
    for (size_t k=0; k<nodes.size(); k++)
    {
        for (size_t j=0; j<nodes.size(); j++)
        {
            for (size_t i=0; i<nodes.size(); i++)
            {
                result.push_back(Nx[i] * Ny[j]  * Nz[k]);
            }
        }
    }
    return result;
}

std::vector<Eigen::Vector3d> DerivativeShapeFunctionsBrickLagrange(const Eigen::Vector3d x, const std::vector<double>& nodes)
{
    std::vector<double> DNx = ShapeFunctions1D::DerivativeShapeFunctionsTrussLagrange(x(0), nodes);
    std::vector<double> DNy = ShapeFunctions1D::DerivativeShapeFunctionsTrussLagrange(x(1), nodes);
    std::vector<double> DNz = ShapeFunctions1D::DerivativeShapeFunctionsTrussLagrange(x(2), nodes);

    std::vector<double> Nx = ShapeFunctions1D::ShapeFunctionsTrussLagrange(x(0), nodes);
    std::vector<double> Ny = ShapeFunctions1D::ShapeFunctionsTrussLagrange(x(1), nodes);
    std::vector<double> Nz = ShapeFunctions1D::ShapeFunctionsTrussLagrange(x(2), nodes);

    std::vector<Eigen::Vector3d> result;
    for (size_t k=0; k<nodes.size(); k++)
    {
        for (size_t j=0; j<nodes.size(); j++)
        {
            for (size_t i=0; i<nodes.size(); i++)
            {
                result.push_back( Eigen::Vector3d({DNx[i]*Ny[j]*Nz[k] , Nx[i]*DNy[j]*Nz[k], Nx[i]*Ny[j]*DNz[k]})   );
            }
        }
    }
    return result;
}

}

} /* namespace NuTo */
