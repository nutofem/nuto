#pragma once

#include <vector>
#include "nuto/mechanics/nodes/CoordinateNode.h"
#include "nuto/mechanics/elements/ElementInterface.h"
#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/mechanics/cell/Matrix.h"

namespace NuTo
{
class CoordinateElementFem : public ElementInterface
{
public:
    CoordinateElementFem(std::vector<CoordinateNode*> nodes, const InterpolationSimple& interpolation)
        : mInterpolation(interpolation)
        , mShape(interpolation.GetShape())
    {
        for (CoordinateNode* node : nodes)
            mNodes.push_back(*node);
        assert(static_cast<int>(mNodes.size()) == interpolation.GetNumNodes());
    }

    CoordinateElementFem(std::initializer_list<std::reference_wrapper<CoordinateNode>> nodes, const InterpolationSimple& interpolation)
        : mNodes(nodes)
        , mInterpolation(interpolation)
        , mShape(interpolation.GetShape())
    {
        assert(static_cast<int>(mNodes.size()) == interpolation.GetNumNodes());
    }

    virtual Eigen::VectorXd ExtractNodeValues(int instance = 0) const override
    {
        // Solve this later, by using type traits
        assert(instance == 0 && "Coordinate nodes can have only 1 instance");

        const int dim = GetDofDimension();
        Eigen::VectorXd nodeValues(GetNumNodes() * dim);
        for (int i = 0; i < GetNumNodes(); ++i)
            nodeValues.segment(dim * i, dim) = GetNode(i).GetCoordinates();
        return nodeValues;
    }

    virtual Eigen::MatrixXd GetNMatrix(NaturalCoords ipCoords) const override
    {
        return Matrix::N(Interpolation().GetShapeFunctions(ipCoords), Interpolation().GetNumNodes(), GetDofDimension());
    }

    virtual Eigen::VectorXd GetShapeFunctions(NaturalCoords ipCoords) const override
    {
        return Interpolation().GetShapeFunctions(ipCoords);
    }

    virtual Eigen::MatrixXd GetDerivativeShapeFunctions(NaturalCoords ipCoords) const override
    {
        return Interpolation().GetDerivativeShapeFunctions(ipCoords);
    }

    virtual int GetDofDimension() const override
    {
        return GetNode(0).GetNumValues();
    }

    Eigen::VectorXi GetDofNumbering() const override
    {
        assert(false && "Coordinate nodes have no dof numbering");
    }

    virtual int GetNumNodes() const override
    {
        return mNodes.size();
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
    const Shape& mShape;
};

} /* NuTo */
