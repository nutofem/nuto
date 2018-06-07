#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/mechanics/interpolation/InterpolationTrussLobatto.h"
#include "nuto/mechanics/interpolation/InterpolationBrickLobatto.h"
#include "nuto/math/shapes/Hexahedron.h"

namespace NuTo
{
Eigen::MatrixXd InterpolationBrickLobatto::LocalCoords(int nodeId, const Eigen::VectorXd& nodes)
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


Eigen::VectorXd InterpolationBrickLobatto::ShapeFunctions(const Eigen::Vector3d x, const Eigen::VectorXd& nodes)
{
    Eigen::VectorXd Nx = InterpolationTrussLobatto::ShapeFunctions(x(0), nodes);
    Eigen::VectorXd Ny = InterpolationTrussLobatto::ShapeFunctions(x(1), nodes);
    Eigen::VectorXd Nz = InterpolationTrussLobatto::ShapeFunctions(x(2), nodes);

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


Eigen::MatrixXd InterpolationBrickLobatto::DerivativeShapeFunctions(const Eigen::Vector3d x,
                                                                    const Eigen::VectorXd& nodes)
{
    Eigen::VectorXd DNx = InterpolationTrussLobatto::DerivativeShapeFunctions(x(0), nodes);
    Eigen::VectorXd DNy = InterpolationTrussLobatto::DerivativeShapeFunctions(x(1), nodes);
    Eigen::VectorXd DNz = InterpolationTrussLobatto::DerivativeShapeFunctions(x(2), nodes);

    Eigen::VectorXd Nx = InterpolationTrussLobatto::ShapeFunctions(x(0), nodes);
    Eigen::VectorXd Ny = InterpolationTrussLobatto::ShapeFunctions(x(1), nodes);
    Eigen::VectorXd Nz = InterpolationTrussLobatto::ShapeFunctions(x(2), nodes);

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

InterpolationBrickLobatto::InterpolationBrickLobatto(int order)
{
    mNodes = InterpolationTrussLobatto::LocalCoords(order);
}

std::unique_ptr<InterpolationSimple> InterpolationBrickLobatto::Clone() const
{
    return std::make_unique<InterpolationBrickLobatto>(*this);
}

Eigen::VectorXd InterpolationBrickLobatto::GetShapeFunctions(const NaturalCoords& naturalIpCoords) const
{
    return ShapeFunctions(naturalIpCoords, mNodes);
}

Eigen::MatrixXd InterpolationBrickLobatto::GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const
{
    return DerivativeShapeFunctions(naturalIpCoords, mNodes);
}

NaturalCoords InterpolationBrickLobatto::GetLocalCoords(int nodeId) const
{
    return LocalCoords(nodeId, mNodes);
}

int InterpolationBrickLobatto::GetNumNodes() const
{
    return mNodes.size() * mNodes.size() * mNodes.size();
}

const Shape& InterpolationBrickLobatto::GetShape() const
{
    return mShape;
}

std::vector<int> InterpolationBrickLobatto::EdgeNodeIds(int edgeIndex) const
{
    throw Exception(__PRETTY_FUNCTION__, "Not implemented");
}

std::unique_ptr<InterpolationSimple> InterpolationBrickLobatto::EdgeInterpolation(int /* edgeIndex*/) const
{
    throw Exception(__PRETTY_FUNCTION__, "Not implemented");
}

std::vector<int> InterpolationBrickLobatto::FaceNodeIds(int /* faceIndex */) const
{
    throw Exception(__PRETTY_FUNCTION__, "Not implemented");
}

std::unique_ptr<InterpolationSimple> InterpolationBrickLobatto::FaceInterpolation(int /* faceIndex*/) const
{
    throw Exception(__PRETTY_FUNCTION__, "Not implemented");
}

} /* NuTo */
