#pragma once
#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/mechanics/interpolation/InterpolationTrussLobatto.h"
#include "nuto/math/shapes/Hexahedron.h"

namespace NuTo
{
class InterpolationBrickLobatto : public InterpolationSimple
{
public:

    static Eigen::MatrixXd NodeCoordinatesBrickLobatto(int nodeId, const Eigen::VectorXd& nodes)
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


    static Eigen::VectorXd ShapeFunctionsBrickLagrange(const Eigen::Vector3d x, const Eigen::VectorXd& nodes)
    {
        Eigen::VectorXd Nx = InterpolationTrussLobatto::ShapeFunctionsTrussLagrange(x(0), nodes);
        Eigen::VectorXd Ny = InterpolationTrussLobatto::ShapeFunctionsTrussLagrange(x(1), nodes);
        Eigen::VectorXd Nz = InterpolationTrussLobatto::ShapeFunctionsTrussLagrange(x(2), nodes);

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


    static Eigen::MatrixXd DerivativeShapeFunctionsBrickLagrange(const Eigen::Vector3d x,
                                                                            const Eigen::VectorXd& nodes)
    {
        Eigen::VectorXd DNx = InterpolationTrussLobatto::DerivativeShapeFunctionsTrussLagrange(x(0), nodes);
        Eigen::VectorXd DNy = InterpolationTrussLobatto::DerivativeShapeFunctionsTrussLagrange(x(1), nodes);
        Eigen::VectorXd DNz = InterpolationTrussLobatto::DerivativeShapeFunctionsTrussLagrange(x(2), nodes);

        Eigen::VectorXd Nx = InterpolationTrussLobatto::ShapeFunctionsTrussLagrange(x(0), nodes);
        Eigen::VectorXd Ny = InterpolationTrussLobatto::ShapeFunctionsTrussLagrange(x(1), nodes);
        Eigen::VectorXd Nz = InterpolationTrussLobatto::ShapeFunctionsTrussLagrange(x(2), nodes);

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

    InterpolationBrickLobatto(int order)
    {
        mNodes = InterpolationTrussLobatto::NodeCoordinatesTrussLobatto(order);
    }

    std::unique_ptr<InterpolationSimple> Clone() const override
    {
        return std::make_unique<InterpolationBrickLobatto>(*this);
    }

    Eigen::VectorXd GetShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        return ShapeFunctionsBrickLagrange(naturalIpCoords, mNodes);
    }

    DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        return DerivativeShapeFunctionsBrickLagrange(naturalIpCoords, mNodes);
    }

    NaturalCoords GetLocalCoords(int nodeId) const override
    {
        return NodeCoordinatesBrickLobatto(nodeId, mNodes);
    }

    int GetNumNodes() const override
    {
        return mNodes.size() * mNodes.size() * mNodes.size();
    }

    const Shape& GetShape() const override
    {
        return mShape;
    }

private:
    Eigen::VectorXd mNodes;
    Hexahedron mShape;
};
} /* NuTo */
