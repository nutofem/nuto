#pragma once

#include <vector>
#include "nuto/mechanics/nodes/DofNode.h"
#include "nuto/mechanics/elements/DofElementInterface.h"
#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/mechanics/cell/Matrix.h"

namespace NuTo
{
class DofElementFem : public DofElementInterface
{
public:
    DofElementFem(std::vector<DofNode*> nodes, const InterpolationSimple& interpolation)
        : mInterpolation(interpolation)
        , mShape(interpolation.GetShape())
    {
        for (DofNode* node : nodes)
            mNodes.push_back(*node);
        assert(static_cast<int>(mNodes.size()) == interpolation.GetNumNodes());
    }

    DofElementFem(std::initializer_list<std::reference_wrapper<DofNode>> nodes,
                  const InterpolationSimple& interpolation)
        : mNodes(nodes)
        , mInterpolation(interpolation)
        , mShape(interpolation.GetShape())
    {
        assert(static_cast<int>(mNodes.size()) == interpolation.GetNumNodes());
    }

    Eigen::VectorXd ExtractNodeValues(int instance = 0) const override
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

    DofNode& GetNode(int i)
    {
        assert(i < static_cast<int>(mNodes.size()));
        return mNodes[i];
    }


    const DofNode& GetNode(int i) const
    {
        assert(i < static_cast<int>(mNodes.size()));
        return mNodes[i];
    }

    const Shape& GetShape() const
    {
        return mShape;
    }

private:
    std::vector<std::reference_wrapper<DofNode>> mNodes;
    std::reference_wrapper<const InterpolationSimple> mInterpolation;
    const Shape& mShape;
};

} /* NuTo */
