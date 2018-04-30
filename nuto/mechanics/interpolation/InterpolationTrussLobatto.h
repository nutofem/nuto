#pragma once
#include <eigen3/Eigen/Core>
#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/math/Legendre.h"
#include "nuto/math/shapes/Line.h"
#include "nuto/base/Exception.h"

namespace NuTo
{
class InterpolationTrussLobatto : public InterpolationSimple
{
public:

    static Eigen::VectorXd NodeCoordinatesTrussLobatto(int order)
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


    static Eigen::VectorXd BarycentricWeights(const Eigen::VectorXd& nodes)
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


    static Eigen::VectorXd ShapeFunctionsTrussLagrange(const double x, const Eigen::VectorXd& nodes)
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


    static Eigen::VectorXd DerivativeShapeFunctionsTrussLagrange(const double x, const Eigen::VectorXd& nodes)
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

    InterpolationTrussLobatto(int order)
    {
        mNodes = NodeCoordinatesTrussLobatto(order);
    }

    std::unique_ptr<InterpolationSimple> Clone() const override
    {
        return std::make_unique<InterpolationTrussLobatto>(*this);
    }

    Eigen::VectorXd GetShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        Eigen::VectorXd result(mNodes.size());
        const Eigen::VectorXd shapes = ShapeFunctionsTrussLagrange(naturalIpCoords[0], mNodes);
        for (int i = 0; i < mNodes.size(); i++)
        {
            result[i] = shapes[i];
        }
        return result;
    }

    DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        Eigen::VectorXd result(mNodes.size());
        const Eigen::VectorXd shapes =
                DerivativeShapeFunctionsTrussLagrange(naturalIpCoords[0], mNodes);
        for (int i = 0; i < mNodes.size(); i++)
        {
            result[i] = shapes[i];
        }
        return result;
    }

    NaturalCoords GetLocalCoords(int nodeId) const override
    {
        Eigen::VectorXd result( 1);
        result[0] = mNodes[nodeId];
        return result;
    }

    int GetNumNodes() const override
    {
        return mNodes.size();
    }

    const Shape& GetShape() const override
    {
        return mShape;
    }

private:
    Eigen::VectorXd mNodes;
    Line mShape;
};
} /* NuTo */
