#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/mechanics/interpolation/InterpolationTrussLobatto.h"
#include "nuto/mechanics/interpolation/InterpolationQuadLobatto.h"
#include "nuto/math/shapes/Quadrilateral.h"

namespace NuTo
{
     Eigen::MatrixXd InterpolationQuadLobatto::LocalCoords(int nodeId, const Eigen::VectorXd& nodes)
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

    Eigen::VectorXd InterpolationQuadLobatto::ShapeFunctions(const Eigen::Vector2d x, const Eigen::VectorXd& nodes)
    {
        Eigen::VectorXd Nx = InterpolationTrussLobatto::ShapeFunctions(x(0), nodes);
        Eigen::VectorXd Ny = InterpolationTrussLobatto::ShapeFunctions(x(1), nodes);

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

    Eigen::MatrixXd InterpolationQuadLobatto::DerivativeShapeFunctions(const Eigen::Vector2d x,
                                                                           const Eigen::VectorXd& nodes)
    {
        Eigen::VectorXd DNx = InterpolationTrussLobatto::DerivativeShapeFunctions(x(0), nodes);
        Eigen::VectorXd DNy = InterpolationTrussLobatto::DerivativeShapeFunctions(x(1), nodes);

        Eigen::VectorXd Nx = InterpolationTrussLobatto::ShapeFunctions(x(0), nodes);
        Eigen::VectorXd Ny = InterpolationTrussLobatto::ShapeFunctions(x(1), nodes);

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

    InterpolationQuadLobatto::InterpolationQuadLobatto(int order)
    {
        mNodes = InterpolationTrussLobatto::LocalCoords(order);
    }

    std::unique_ptr<InterpolationSimple> InterpolationQuadLobatto::Clone() const
    {
        return std::make_unique<InterpolationQuadLobatto>(*this);
    }

    Eigen::VectorXd InterpolationQuadLobatto::GetShapeFunctions(const NaturalCoords& naturalIpCoords) const
    {
        return ShapeFunctions(naturalIpCoords, mNodes);
    }

    DerivativeShapeFunctionsNatural InterpolationQuadLobatto::GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const
    {
        return DerivativeShapeFunctions(naturalIpCoords, mNodes);
    }

    NaturalCoords InterpolationQuadLobatto::GetLocalCoords(int nodeId) const
    {
        return LocalCoords(nodeId, mNodes);
    }

    int InterpolationQuadLobatto::GetNumNodes() const
    {
        return mNodes.size() * mNodes.size();
    }

    const Shape& InterpolationQuadLobatto::GetShape() const
    {
        return mShape;
    }

} /* NuTo */
