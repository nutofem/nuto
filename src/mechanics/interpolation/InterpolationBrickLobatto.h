#pragma once
#include "mechanics/interpolation/InterpolationSimple.h"
#include "mechanics/elements/SpectralShapeFunctions.h"

namespace NuTo
{
class InterpolationBrickLobatto : public InterpolationSimple
{
public:
    InterpolationBrickLobatto(int order)
    {
        mNodes = ShapeFunctions1D::NodeCoordinatesTrussLobatto(order);
    }

    std::unique_ptr<InterpolationSimple> Clone() const override
    {
        return std::make_unique<InterpolationBrickLobatto>(*this);
    }

    ShapeFunctions GetShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        return ShapeFunctions3D::ShapeFunctionsBrickLagrange(naturalIpCoords, mNodes);
    }

    DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        return ShapeFunctions3D::DerivativeShapeFunctionsBrickLagrange(naturalIpCoords, mNodes);
    }

    NaturalCoords GetLocalCoords(int nodeId) const override
    {
        return ShapeFunctions3D::NodeCoordinatesBrickLobatto(nodeId, mNodes);
    }

    int GetNumNodes() const override
    {
        return mNodes.size() * mNodes.size() * mNodes.size();
    }

private:
    Eigen::VectorXd mNodes;
};
} /* NuTo */
