#include "SpectralShapeFunctions.h"
#include "math/Legendre.h"
#include "base/Exception.h"

using namespace NuTo;

Eigen::VectorXd ShapeFunctions1D::NodeCoordinatesTrussLobatto(int order)
{
    if (order < 1)
        throw NuTo::Exception(__PRETTY_FUNCTION__, "Order too low. Must be 1 or higher");

    std::vector<double> points = NuTo::Math::Polynomial::LegendreDerivRoots(order);

    Eigen::VectorXd result(points.size() + 2);
    result[0] = -1.;
    for (size_t i = 0; i < points.size(); i++)
    {
        result[i + 1] = points[i];
    }
    result[result.rows() - 1] = 1.;
    return result;
}

Eigen::VectorXd ShapeFunctions1D::BarycentricWeights(const Eigen::VectorXd& nodes)
{
    Eigen::VectorXd w = Eigen::VectorXd::Ones(nodes.rows());
    for (int j = 0; j < nodes.rows(); j++)
        for (int i = 0; i < nodes.rows(); i++)
        {
            if (i != j)
            {
                w[j] /= (nodes[j] - nodes[i]);
            }
        }
    return w;
}


Eigen::VectorXd ShapeFunctions1D::ShapeFunctionsTrussLagrange(const double x, const Eigen::VectorXd& nodes)
{
    Eigen::VectorXd result = Eigen::VectorXd::Zero(nodes.rows());
    Eigen::VectorXd w = BarycentricWeights(nodes);
    // Check if x is a node (or near)
    bool xMatchesNode = false;
    for (int j = 0; j < nodes.rows(); j++)
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

    for (int j = 0; j < nodes.rows(); j++)
    {
        double tmp = w[j] / (x - nodes[j]);
        result[j] = tmp;
        sum += tmp;
    }

    for (int j = 0; j < nodes.rows(); j++)
    {
        result[j] /= sum;
    }
    return result;
}


Eigen::VectorXd ShapeFunctions1D::DerivativeShapeFunctionsTrussLagrange(const double x, const Eigen::VectorXd& nodes)
{
    Eigen::VectorXd result = Eigen::VectorXd::Zero(nodes.rows());
    Eigen::VectorXd w = BarycentricWeights(nodes);
    for (int j = 0; j < nodes.rows(); j++)
    {
        for (int k = 0; k < nodes.rows(); k++)
        {
            if (k != j)
            {
                double tmp = 1.;
                for (int i = 0; i < nodes.rows(); i++)
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

Eigen::MatrixXd ShapeFunctions2D::NodeCoordinatesQuadLobatto(int nodeId, const Eigen::VectorXd& nodes)
{
    const int d = nodes.rows();

    assert(nodeId >= 0);
    assert(nodeId < nodes.rows() * nodes.rows());

    int i = nodeId % d;
    int j = nodeId / d;

    double cX = nodes[i];
    double cY = nodes[j];

    return Eigen::Vector2d({cX, cY});
}

Eigen::VectorXd ShapeFunctions2D::ShapeFunctionsQuadLagrange(const Eigen::Vector2d x, const Eigen::VectorXd& nodes)
{
    Eigen::VectorXd Nx = ShapeFunctions1D::ShapeFunctionsTrussLagrange(x(0), nodes);
    Eigen::VectorXd Ny = ShapeFunctions1D::ShapeFunctionsTrussLagrange(x(1), nodes);

    Eigen::VectorXd result(nodes.rows() * nodes.rows());
    int count = 0;
    for (int j = 0; j < nodes.rows(); j++)
    {
        for (int i = 0; i < nodes.rows(); i++)
        {
            result[count] = Nx[i] * Ny[j];
            count++;
        }
    }
    return result;
}

Eigen::MatrixXd ShapeFunctions2D::DerivativeShapeFunctionsQuadLagrange(const Eigen::Vector2d x,
                                                                       const Eigen::VectorXd& nodes)
{
    Eigen::VectorXd DNx = ShapeFunctions1D::DerivativeShapeFunctionsTrussLagrange(x(0), nodes);
    Eigen::VectorXd DNy = ShapeFunctions1D::DerivativeShapeFunctionsTrussLagrange(x(1), nodes);

    Eigen::VectorXd Nx = ShapeFunctions1D::ShapeFunctionsTrussLagrange(x(0), nodes);
    Eigen::VectorXd Ny = ShapeFunctions1D::ShapeFunctionsTrussLagrange(x(1), nodes);

    Eigen::MatrixXd result = Eigen::MatrixXd::Zero(nodes.rows() * nodes.rows(), 2);
    int count = 0;
    for (int j = 0; j < nodes.rows(); j++)
    {
        for (int i = 0; i < nodes.rows(); i++)
        {
            result(count, 0) = DNx[i] * Ny[j];
            result(count, 1) = Nx[i] * DNy[j];
            count++;
        }
    }
    return result;
}


Eigen::MatrixXd ShapeFunctions3D::NodeCoordinatesBrickLobatto(int nodeId, const Eigen::VectorXd& nodes)
{
    const int d = nodes.rows();

    assert(nodeId >= 0);
    assert(nodeId < nodes.rows() * nodes.rows() * nodes.rows());

    int i = nodeId % d;
    int j = nodeId % (d * d) / d;
    int k = nodeId / (d * d);

    double cX = nodes[i];
    double cY = nodes[j];
    double cZ = nodes[k];

    return Eigen::Vector3d({cX, cY, cZ});
}


Eigen::VectorXd ShapeFunctions3D::ShapeFunctionsBrickLagrange(const Eigen::Vector3d x, const Eigen::VectorXd& nodes)
{
    Eigen::VectorXd Nx = ShapeFunctions1D::ShapeFunctionsTrussLagrange(x(0), nodes);
    Eigen::VectorXd Ny = ShapeFunctions1D::ShapeFunctionsTrussLagrange(x(1), nodes);
    Eigen::VectorXd Nz = ShapeFunctions1D::ShapeFunctionsTrussLagrange(x(2), nodes);

    Eigen::VectorXd result = Eigen::VectorXd::Zero(nodes.rows() * nodes.rows() * nodes.rows());
    int count = 0;
    for (int k = 0; k < nodes.rows(); k++)
    {
        for (int j = 0; j < nodes.rows(); j++)
        {
            for (int i = 0; i < nodes.rows(); i++)
            {
                result[count] = Nx[i] * Ny[j] * Nz[k];
                count++;
            }
        }
    }
    return result;
}


Eigen::MatrixXd ShapeFunctions3D::DerivativeShapeFunctionsBrickLagrange(const Eigen::Vector3d x,
                                                                        const Eigen::VectorXd& nodes)
{
    Eigen::VectorXd DNx = ShapeFunctions1D::DerivativeShapeFunctionsTrussLagrange(x(0), nodes);
    Eigen::VectorXd DNy = ShapeFunctions1D::DerivativeShapeFunctionsTrussLagrange(x(1), nodes);
    Eigen::VectorXd DNz = ShapeFunctions1D::DerivativeShapeFunctionsTrussLagrange(x(2), nodes);

    Eigen::VectorXd Nx = ShapeFunctions1D::ShapeFunctionsTrussLagrange(x(0), nodes);
    Eigen::VectorXd Ny = ShapeFunctions1D::ShapeFunctionsTrussLagrange(x(1), nodes);
    Eigen::VectorXd Nz = ShapeFunctions1D::ShapeFunctionsTrussLagrange(x(2), nodes);

    Eigen::MatrixXd result = Eigen::MatrixXd::Zero(nodes.rows() * nodes.rows() * nodes.rows(), 3);
    int count = 0;
    for (int k = 0; k < nodes.rows(); k++)
    {
        for (int j = 0; j < nodes.rows(); j++)
        {
            for (int i = 0; i < nodes.rows(); i++)
            {
                result(count, 0) = DNx[i] * Ny[j] * Nz[k];
                result(count, 1) = Nx[i] * DNy[j] * Nz[k];
                result(count, 2) = Nx[i] * Ny[j] * DNz[k];
                count++;
            }
        }
    }
    return result;
}
