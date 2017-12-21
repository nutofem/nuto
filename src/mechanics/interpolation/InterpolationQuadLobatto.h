#pragma once
#include "mechanics/interpolation/InterpolationSimple.h"
#include "mechanics/elements/SpectralShapeFunctions.h"

namespace NuTo
{
class InterpolationQuadLobatto : public InterpolationSimple
{
public:
    InterpolationQuadLobatto(int dofDimension, int order)
        : mDofDimension(dofDimension)
    {
        mNodes = ShapeFunctions1D::NodeCoordinatesTrussLobatto(order);
    }

    std::unique_ptr<InterpolationSimple> Clone() const override
    {
        return std::make_unique<InterpolationQuadLobatto>(*this);
    }

    ShapeFunctions GetShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        return ShapeFunctions2D::ShapeFunctionsQuadLagrange(naturalIpCoords, mNodes);
    }

    DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        return ShapeFunctions2D::DerivativeShapeFunctionsQuadLagrange(naturalIpCoords, mNodes);
    }

    NaturalCoords GetLocalCoords(int nodeId) const override
    {
        const int d = mNodes.size();

        assert(nodeId >= 0);
        assert(nodeId < GetNumNodes());

        int i = nodeId % d;
        int j = nodeId / d;

        double cX = mNodes[i];
        double cY = mNodes[j];

        return Eigen::Vector2d({cX, cY});
    }

    int GetNumNodes() const override
    {
        return mNodes.size() * mNodes.size();
    }

    int GetDofDimension() const override
    {
        return mDofDimension;
    }

private:
    int mDofDimension;
    std::vector<double> mNodes;
};
} /* NuTo */
