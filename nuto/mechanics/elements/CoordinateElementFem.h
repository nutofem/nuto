#pragma once

#include <vector>
#include "nuto/base/Exception.h"
#include "nuto/mechanics/nodes/CoordinateNode.h"
#include "nuto/mechanics/elements/CoordinateElementInterface.h"
#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/mechanics/cell/Matrix.h"

namespace NuTo
{
class CoordinateElementFem : public CoordinateElementInterface
{
public:
    CoordinateElementFem(std::vector<CoordinateNode*> nodes, const InterpolationSimple& interpolation)
        : mInterpolation(interpolation)
        , mShape(interpolation.GetShape())
    {
        for (CoordinateNode* node : nodes)
            mNodes.push_back(*node);
        assert(static_cast<int>(mNodes.size()) == interpolation.GetNumNodes());
        SetId(-1);
    }

    CoordinateElementFem(std::initializer_list<std::reference_wrapper<CoordinateNode>> nodes,
                         const InterpolationSimple& interpolation)
        : mNodes(nodes)
        , mInterpolation(interpolation)
        , mShape(interpolation.GetShape())
    {
        assert(static_cast<int>(mNodes.size()) == interpolation.GetNumNodes());
        SetId(-1);
    }

    Eigen::VectorXd ExtractCoordinates() const override
    {
        const int dim = GetDofDimension();
        Eigen::VectorXd nodeValues(GetNumNodes() * dim);
        for (int i = 0; i < GetNumNodes(); ++i)
            nodeValues.segment(dim * i, dim) = GetNode(i).GetCoordinates();
        return nodeValues;
    }

    Eigen::MatrixXd GetNMatrix(NaturalCoords ipCoords) const override
    {
        return Matrix::N(Interpolation().GetShapeFunctions(ipCoords), Interpolation().GetNumNodes(), GetDofDimension());
    }

    Eigen::VectorXd GetShapeFunctions(NaturalCoords ipCoords) const override
    {
        return Interpolation().GetShapeFunctions(ipCoords);
    }

    Eigen::MatrixXd GetDerivativeShapeFunctions(NaturalCoords ipCoords) const override
    {
        return Interpolation().GetDerivativeShapeFunctions(ipCoords);
    }

    int GetDofDimension() const override
    {
        return GetNode(0).GetNumValues();
    }

    int GetNumNodes() const override
    {
        return mNodes.size();
    }

    int GetId() const
    {
        return id;
    }

    void SetId(int i)
    {
        id = i;
    }

    const InterpolationSimple& Interpolation() const
    {
        return mInterpolation;
    }

    CoordinateNode& GetNode(int i)
    {
        assert(i < static_cast<int>(mNodes.size()));
        return mNodes[i];
    }


    const CoordinateNode& GetNode(int i) const
    {
        assert(i < static_cast<int>(mNodes.size()));
        return mNodes[i];
    }

    const Shape& GetShape() const
    {
        return mShape;
    }

private:
    std::vector<std::reference_wrapper<CoordinateNode>> mNodes;
    std::reference_wrapper<const InterpolationSimple> mInterpolation;
    int id;
    const Shape& mShape;
};

} /* NuTo */
