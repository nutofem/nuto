#pragma once
#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/mechanics/interpolation/InterpolationTrussLobatto.h"
#include "nuto/math/shapes/Quadrilateral.h"

namespace NuTo
{
class InterpolationQuadLobatto : public InterpolationSimple
{
public:

    static Eigen::MatrixXd NodeCoordinatesQuadLobatto(int nodeId, const Eigen::VectorXd& nodes)
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

    static Eigen::VectorXd ShapeFunctionsQuadLagrange(const Eigen::Vector2d x, const Eigen::VectorXd& nodes)
    {
        Eigen::VectorXd Nx = InterpolationTrussLobatto::ShapeFunctionsTrussLagrange(x(0), nodes);
        Eigen::VectorXd Ny = InterpolationTrussLobatto::ShapeFunctionsTrussLagrange(x(1), nodes);

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

    static Eigen::MatrixXd DerivativeShapeFunctionsQuadLagrange(const Eigen::Vector2d x,
                                                                           const Eigen::VectorXd& nodes)
    {
        Eigen::VectorXd DNx = InterpolationTrussLobatto::DerivativeShapeFunctionsTrussLagrange(x(0), nodes);
        Eigen::VectorXd DNy = InterpolationTrussLobatto::DerivativeShapeFunctionsTrussLagrange(x(1), nodes);

        Eigen::VectorXd Nx = InterpolationTrussLobatto::ShapeFunctionsTrussLagrange(x(0), nodes);
        Eigen::VectorXd Ny = InterpolationTrussLobatto::ShapeFunctionsTrussLagrange(x(1), nodes);

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

    InterpolationQuadLobatto(int order)
    {
        mNodes = InterpolationTrussLobatto::NodeCoordinatesTrussLobatto(order);
    }

    std::unique_ptr<InterpolationSimple> Clone() const override
    {
        return std::make_unique<InterpolationQuadLobatto>(*this);
    }

    Eigen::VectorXd GetShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        return ShapeFunctionsQuadLagrange(naturalIpCoords, mNodes);
    }

    DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        return DerivativeShapeFunctionsQuadLagrange(naturalIpCoords, mNodes);
    }

    NaturalCoords GetLocalCoords(int nodeId) const override
    {
        return NodeCoordinatesQuadLobatto(nodeId, mNodes);
    }

    int GetNumNodes() const override
    {
        return mNodes.size() * mNodes.size();
    }

    const Shape& GetShape() const override
    {
        return mShape;
    }

private:
    Eigen::VectorXd mNodes;
    Quadrilateral mShape;
};
} /* NuTo */
