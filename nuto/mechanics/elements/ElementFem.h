#pragma once

#include <vector>
#include "nuto/mechanics/nodes/NodeCoordinates.h"
#include "nuto/mechanics/nodes/NodeSimple.h"
#include "nuto/mechanics/elements/ElementInterface.h"
#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/mechanics/cell/Matrix.h"

namespace NuTo
{
template <typename TNode>
class ElementFem;

typedef ElementFem<NodeCoordinates> CoordinateElementFem;
typedef ElementFem<NodeSimple> DofElementFem;


template <typename TNode>
class ElementFem : public ElementInterface
{
public:
    ElementFem(std::vector<TNode*> nodes, const InterpolationSimple& interpolation)
        : mInterpolation(interpolation)
        , mShape(interpolation.GetShape())
    {
        for (TNode* node : nodes)
            mNodes.push_back(*node);
        assert(static_cast<int>(mNodes.size()) == interpolation.GetNumNodes());
    }

    ElementFem(std::initializer_list<std::reference_wrapper<TNode>> nodes, const InterpolationSimple& interpolation)
        : mNodes(nodes)
        , mInterpolation(interpolation)
        , mShape(interpolation.GetShape())
    {
        assert(static_cast<int>(mNodes.size()) == interpolation.GetNumNodes());
    }

    virtual NodeValues ExtractNodeValues(int instance = 0) const override
    {
        const int dim = GetDofDimension();
        Eigen::VectorXd nodeValues(GetNumNodes() * dim);
        for (int i = 0; i < GetNumNodes(); ++i)
            nodeValues.segment(dim * i, dim) = GetNode(i).GetValues(instance);
        return nodeValues;
    }

    NMatrix GetNMatrix(NaturalCoords ipCoords) const override
    {
        return Matrix::N(Interpolation().GetShapeFunctions(ipCoords), Interpolation().GetNumNodes(), GetDofDimension());
    }

    ShapeFunctions GetShapeFunctions(NaturalCoords ipCoords) const override
    {
        return Interpolation().GetShapeFunctions(ipCoords);
    }

    DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(NaturalCoords ipCoords) const override
    {
        return Interpolation().GetDerivativeShapeFunctions(ipCoords);
    }

    int GetDofDimension() const override
    {
        return GetNode(0).GetNumValues();
    }

    Eigen::VectorXi GetDofNumbering() const override
    {
        Eigen::VectorXi dofNumbering(GetNumNodes() * GetDofDimension());
        int i = 0;
        for (int iNode = 0; iNode < GetNumNodes(); ++iNode)
        {
            const auto& node = GetNode(iNode);
            for (int iDof = 0; iDof < GetDofDimension(); ++iDof)
            {
                dofNumbering[i++] = node.GetDofNumber(iDof);
            }
        }
        return dofNumbering;
    }

    int GetNumNodes() const override
    {
        return mNodes.size();
    }

    const InterpolationSimple& Interpolation() const
    {
        return mInterpolation;
    }

    TNode& GetNode(int i)
    {
        assert(i < static_cast<int>(mNodes.size()));
        return mNodes[i];
    }


    const TNode& GetNode(int i) const
    {
        assert(i < static_cast<int>(mNodes.size()));
        return mNodes[i];
    }

    const Shape& GetShape() const
    {
        return mShape;
    }

private:
    std::vector<std::reference_wrapper<TNode>> mNodes;
    std::reference_wrapper<const InterpolationSimple> mInterpolation;
    const Shape& mShape;
};
} /* NuTo */
