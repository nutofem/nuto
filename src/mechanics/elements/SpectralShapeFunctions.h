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
Eigen::VectorXd NodeCoordinatesTrussLobatto(int order)
{
    if (order < 1)
        throw NuTo::Exception(__PRETTY_FUNCTION__, "Order too low. Must be 1 or higher");

    std::vector<double> points = NuTo::Math::Polynomial::LegendreDerivRoots(order);

    Eigen::VectorXd result(points.size()+2);
    result[0] = -1.;
    for (size_t i = 0; i<points.size();i++)
    {
        result[i+1] = points[i];
    }
    result[result.size()-1] = 1.;
    return result;
}

//! @brief Compute barycentric weights from node coordinates
//! @param nodes local node coordinates
Eigen::VectorXd BarycentricWeights(const Eigen::VectorXd& nodes)
{
    Eigen::VectorXd w = Eigen::VectorXd::Ones(nodes.size());
    for (int j = 0; j < nodes.size(); j++)
        for (int i = 0; i < nodes.size(); i++)
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
Eigen::VectorXd ShapeFunctionsTrussLagrange(const double x, const Eigen::VectorXd& nodes)
{
    Eigen::VectorXd result = Eigen::VectorXd::Zero(nodes.size());
    Eigen::VectorXd w = BarycentricWeights(nodes);
    // Check if x is a node (or near)
    bool xMatchesNode = false;
    for (int j = 0; j < nodes.size(); j++)
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

    for (int j = 0; j < nodes.size(); j++)
    {
        double tmp = w[j]/(x-nodes[j]);
        result[j] = tmp;
        sum += tmp;
    }

    for (int j = 0; j < nodes.size(); j++)
    {
        result[j] /= sum;
    }
    return result;
}

//! @brief Value of 1D Lagrange type shape function derivative with given nodes
//! @param x local coordinate where shapes are evaluated
//! @param nodes local node coordinates
//! @return shape derivatives evaluated at x.
Eigen::VectorXd DerivativeShapeFunctionsTrussLagrange(const double x, const Eigen::VectorXd& nodes)
{
    Eigen::VectorXd result = Eigen::VectorXd::Zero(nodes.size());
    Eigen::VectorXd w = BarycentricWeights(nodes);
    for (int j = 0; j < nodes.size(); j++)
    {
        for (int k = 0; k < nodes.size(); k++)
        {
            if (k != j)
            {
                double tmp = 1.;
                for (int i = 0; i < nodes.size(); i++)
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
Eigen::MatrixXd NodeCoordinatesQuadLobatto(int order)
{
    Eigen::VectorXd nodes = ShapeFunctions1D::NodeCoordinatesTrussLobatto(order);
    Eigen::MatrixXd result(nodes.size() * nodes.size(),2);
    int count = 0;
    for (int i=0; i<nodes.size(); i++)
    {
        for (int j=0; j<nodes.size(); j++)
        {
            result(count,0) = nodes[i];
            result(count,1) = nodes[j];
            count++;
        }
    }
    return result;
}

//! @brief Value of Lagrange type shape functions with given 1D nodes (tensor product)
//! @param x local coordinate where shapes are evaluated
//! @param nodes local node coordinates
//! @return shapes evaluated at x. Fulfill interpolation condition fi(xj) = delta_ij
Eigen::VectorXd ShapeFunctionsQuadLagrange(const Eigen::Vector2d x, const Eigen::VectorXd& nodes)
{
    Eigen::VectorXd Nx = ShapeFunctions1D::ShapeFunctionsTrussLagrange(x(0), nodes);
    Eigen::VectorXd Ny = ShapeFunctions1D::ShapeFunctionsTrussLagrange(x(1), nodes);

    Eigen::VectorXd result(nodes.size() * nodes.size());
    int count = 0;
    for (int j=0; j<nodes.size(); j++)
    {
        for (int i=0; i<nodes.size(); i++)
        {
            result[count] = Nx[i] * Ny[j];
            count++;
        }
    }
    return result;
}

Eigen::MatrixXd DerivativeShapeFunctionsQuadLagrange(const Eigen::Vector2d x, const Eigen::VectorXd& nodes)
{
    Eigen::VectorXd DNx = ShapeFunctions1D::DerivativeShapeFunctionsTrussLagrange(x(0), nodes);
    Eigen::VectorXd DNy = ShapeFunctions1D::DerivativeShapeFunctionsTrussLagrange(x(1), nodes);

    Eigen::VectorXd Nx = ShapeFunctions1D::ShapeFunctionsTrussLagrange(x(0), nodes);
    Eigen::VectorXd Ny = ShapeFunctions1D::ShapeFunctionsTrussLagrange(x(1), nodes);

    Eigen::MatrixXd result = Eigen::MatrixXd::Zero(nodes.size() * nodes.size(),2);
    int count = 0;
    for (int j=0; j<nodes.size(); j++)
    {
        for (int i=0; i<nodes.size(); i++)
        {
            result(count,0) = DNx[i]*Ny[j];
            result(count,1) =  Nx[i]*DNy[j];
            count++;
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
Eigen::MatrixXd NodeCoordinatesBrickLobatto(int order)
{
    Eigen::VectorXd nodes = ShapeFunctions1D::NodeCoordinatesTrussLobatto(order);
    Eigen::MatrixXd result(nodes.size() * nodes.size() * nodes.size(),3);
    int count = 0;
    for (int i=0; i<nodes.size(); i++)
    {
        for (int j=0; j<nodes.size(); j++)
        {
            for (int k=0; k<nodes.size(); k++)
            {
                result(count,0) = nodes[i];
                result(count,1) = nodes[j];
                result(count,2) = nodes[k];
                count++;
            }
        }
    }
    return result;
}

//! @brief Value of Lagrange type shape functions with given 1D nodes (tensor product)
//! @param x local coordinate where shapes are evaluated
//! @param nodes local node coordinates
//! @return shapes evaluated at x. Fulfill interpolation condition fi(xj) = delta_ij
Eigen::VectorXd ShapeFunctionsBrickLagrange(const Eigen::Vector3d x, const Eigen::VectorXd& nodes)
{
    Eigen::VectorXd Nx = ShapeFunctions1D::ShapeFunctionsTrussLagrange(x(0), nodes);
    Eigen::VectorXd Ny = ShapeFunctions1D::ShapeFunctionsTrussLagrange(x(1), nodes);
    Eigen::VectorXd Nz = ShapeFunctions1D::ShapeFunctionsTrussLagrange(x(2), nodes);

    Eigen::VectorXd result = Eigen::VectorXd::Zero(nodes.size() * nodes.size() * nodes.size());
    int count = 0;
    for (int k=0; k<nodes.size(); k++)
    {
        for (int j=0; j<nodes.size(); j++)
        {
            for (int i=0; i<nodes.size(); i++)
            {
                result[count] = Nx[i] * Ny[j]  * Nz[k];
                count++;
            }
        }
    }
    return result;
}

Eigen::MatrixXd DerivativeShapeFunctionsBrickLagrange(const Eigen::Vector3d x, const Eigen::VectorXd& nodes)
{
    Eigen::VectorXd DNx = ShapeFunctions1D::DerivativeShapeFunctionsTrussLagrange(x(0), nodes);
    Eigen::VectorXd DNy = ShapeFunctions1D::DerivativeShapeFunctionsTrussLagrange(x(1), nodes);
    Eigen::VectorXd DNz = ShapeFunctions1D::DerivativeShapeFunctionsTrussLagrange(x(2), nodes);

    Eigen::VectorXd Nx = ShapeFunctions1D::ShapeFunctionsTrussLagrange(x(0), nodes);
    Eigen::VectorXd Ny = ShapeFunctions1D::ShapeFunctionsTrussLagrange(x(1), nodes);
    Eigen::VectorXd Nz = ShapeFunctions1D::ShapeFunctionsTrussLagrange(x(2), nodes);

    Eigen::MatrixXd result = Eigen::MatrixXd::Zero(nodes.size() * nodes.size() * nodes.size(), 3);
    int count = 0;
    for (int k=0; k<nodes.size(); k++)
    {
        for (int j=0; j<nodes.size(); j++)
        {
            for (int i=0; i<nodes.size(); i++)
            {
                result(count,0) = DNx[i]* Ny[j]* Nz[k];
                result(count,1) =  Nx[i]*DNy[j]* Nz[k];
                result(count,2) =  Nx[i]* Ny[j]*DNz[k];
                count++;
            }
        }
    }
    return result;
}

}

} /* namespace NuTo */
