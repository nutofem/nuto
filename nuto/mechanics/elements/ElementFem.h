#pragma once

#include <vector>
#include "nuto/mechanics/nodes/CoordinateNode.h"
#include "nuto/mechanics/nodes/DofNode.h"
#include "nuto/mechanics/elements/ElementInterface.h"
#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/mechanics/cell/Matrix.h"

namespace NuTo
{
template <typename TNode>
class ElementFem;

typedef ElementFem<CoordinateNode> CoordinateElementFem;
typedef ElementFem<DofNode> DofElementFem;


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

    virtual Eigen::VectorXd ExtractNodeValues(int instance = 0) const override
    {
        const int dim = GetDofDimension();
        Eigen::VectorXd nodeValues(GetNumNodes() * dim);
        for (int i = 0; i < GetNumNodes(); ++i)
            nodeValues.segment(dim * i, dim) = GetNode(i).GetValues(instance);
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

    Eigen::VectorXi GetDofNumbering() const override;

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

template <>
inline Eigen::VectorXd ElementFem<DofNode>::ExtractNodeValues(int instance) const
{
    const int dim = GetDofDimension();
    Eigen::VectorXd nodeValues(GetNumNodes() * dim);
    for (int i = 0; i < GetNumNodes(); ++i)
        nodeValues.segment(dim * i, dim) = GetNode(i).GetValues(instance);
    return nodeValues;
}

template <>
inline Eigen::VectorXd ElementFem<CoordinateNode>::ExtractNodeValues(int instance) const
{
    // Solve this later, by using type traits
    assert(instance == 0 && "Coordinate nodes can have only 1 instance");


    const int dim = GetDofDimension();
    Eigen::VectorXd nodeValues(GetNumNodes() * dim);
    for (int i = 0; i < GetNumNodes(); ++i)
        nodeValues.segment(dim * i, dim) = GetNode(i).GetCoordinates();
    return nodeValues;
}

template <>
inline Eigen::VectorXi ElementFem<CoordinateNode>::GetDofNumbering() const
{
    assert(false && "Coordinate nodes have no dof numbering");
}

template <>
inline Eigen::VectorXi ElementFem<DofNode>::GetDofNumbering() const
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
} /* NuTo */
