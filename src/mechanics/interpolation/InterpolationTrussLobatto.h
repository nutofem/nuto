#pragma once
#include <eigen3/Eigen/Core>
#include "mechanics/interpolation/InterpolationSimple.h"
#include "mechanics/elements/SpectralShapeFunctions.h"

namespace NuTo
{
class InterpolationTrussLobatto : public InterpolationSimple
{
public:
    InterpolationTrussLobatto(int dofDimension, int order)
        : mDofDimension(dofDimension)
    {
        mNodes = ShapeFunctions1D::NodeCoordinatesTrussLobatto(order);
    }

    std::unique_ptr<InterpolationSimple> Clone() const override
    {
        return std::make_unique<InterpolationTrussLobatto>(*this);
    }

    ShapeFunctions GetShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        Eigen::VectorXd result(mNodes.size());
        const std::vector<double> shapes = ShapeFunctions1D::ShapeFunctionsTrussLagrange(naturalIpCoords[0], mNodes);
        for (size_t i = 0; i < mNodes.size(); i++)
        {
            result[i] = shapes[i];
        }
        return result;
    }

    DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        Eigen::VectorXd result(mNodes.size());
        const std::vector<double> shapes =
                ShapeFunctions1D::DerivativeShapeFunctionsTrussLagrange(naturalIpCoords[0], mNodes);
        for (size_t i = 0; i < mNodes.size(); i++)
        {
            result[i] = shapes[i];
        }
        return result;
    }

    NaturalCoords GetLocalCoords(int nodeId) const override
    {
        Eigen::VectorXd result(1);
        result[0] = mNodes[nodeId];
        return result;
    }

    int GetNumNodes() const override
    {
        return mNodes.size();
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
