#include <eigen3/Eigen/Core>
#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/mechanics/interpolation/InterpolationTrussLobatto.h"
#include "nuto/math/Legendre.h"
#include "nuto/math/shapes/Line.h"
#include "nuto/base/Exception.h"

namespace NuTo
{
Eigen::VectorXd InterpolationTrussLobatto::LocalCoords(int order)
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


Eigen::VectorXd InterpolationTrussLobatto::BarycentricWeights(const Eigen::VectorXd& nodes)
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


Eigen::VectorXd InterpolationTrussLobatto::ShapeFunctions(const double x, const Eigen::VectorXd& nodes)
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


Eigen::VectorXd InterpolationTrussLobatto::DerivativeShapeFunctions(const double x, const Eigen::VectorXd& nodes)
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

InterpolationTrussLobatto::InterpolationTrussLobatto(int order)
{
    mNodes = LocalCoords(order);
}

std::unique_ptr<InterpolationSimple> InterpolationTrussLobatto::Clone() const
{
    return std::make_unique<InterpolationTrussLobatto>(*this);
}

Eigen::VectorXd InterpolationTrussLobatto::GetShapeFunctions(const NaturalCoords& naturalIpCoords) const
{
    Eigen::VectorXd result(mNodes.size());
    const Eigen::VectorXd shapes = ShapeFunctions(naturalIpCoords[0], mNodes);
    for (int i = 0; i < mNodes.size(); i++)
    {
        result[i] = shapes[i];
    }
    return result;
}

Eigen::MatrixXd InterpolationTrussLobatto::GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const
{
    Eigen::VectorXd result(mNodes.size());
    const Eigen::VectorXd shapes = DerivativeShapeFunctions(naturalIpCoords[0], mNodes);
    for (int i = 0; i < mNodes.size(); i++)
    {
        result[i] = shapes[i];
    }
    return result;
}

NaturalCoords InterpolationTrussLobatto::GetLocalCoords(int nodeId) const
{
    Eigen::VectorXd result(1);
    result[0] = mNodes[nodeId];
    return result;
}

int InterpolationTrussLobatto::GetNumNodes() const
{
    return mNodes.size();
}

const Shape& InterpolationTrussLobatto::GetShape() const
{
    return mShape;
}

std::vector<int> InterpolationTrussLobatto::EdgeNodeIds(int edgeIndex) const
{
    throw Exception(__PRETTY_FUNCTION__, "Not implemented");
}

std::unique_ptr<InterpolationSimple> InterpolationTrussLobatto::EdgeInterpolation(int /* edgeIndex*/) const
{
    throw Exception(__PRETTY_FUNCTION__, "Not implemented");
}

std::vector<int> InterpolationTrussLobatto::FaceNodeIds(int /* faceIndex */) const
{
    throw Exception(__PRETTY_FUNCTION__, "Not implemented");
}

std::unique_ptr<InterpolationSimple> InterpolationTrussLobatto::FaceInterpolation(int /* faceIndex*/) const
{
    throw Exception(__PRETTY_FUNCTION__, "Not implemented");
}

} /* NuTo */
