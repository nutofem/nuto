#pragma once
#include "mechanics/interpolation/InterpolationSimple.h"
#include "mechanics/elements/SpectralShapeFunctions.h"

namespace NuTo
{
class InterpolationBrickLobatto : public InterpolationSimple
{
public:
    InterpolationBrickLobatto(int dofDimension, int order)
        : mDofDimension(dofDimension)
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
        const int d = mNodes.size();

        assert(nodeId >= 0);
        assert(nodeId < GetNumNodes());

        int i = nodeId % d;
        int j = nodeId % (d*d) / d;
        int k = nodeId / (d*d);

        double cX = mNodes[i];
        double cY = mNodes[j];
        double cZ = mNodes[k];

        return Eigen::Vector3d({cX, cY,cZ});
    }

    int GetNumNodes() const override
    {
        return mNodes.size() * mNodes.size() * mNodes.size();
    }

    int GetDofDimension() const override
    {
        return mDofDimension;
    }

private:
    int mDofDimension;
    Eigen::VectorXd mNodes;
};
} /* NuTo */
