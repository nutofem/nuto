#pragma once

#include <vector>
#include "mechanics/nodes/NodeSimple.h"
#include "mechanics/elements/ElementInterface.h"
#include "mechanics/interpolation/InterpolationSimple.h"
#include "mechanics/cell/Matrix.h"

namespace NuTo
{

class ElementFem : public ElementInterface
{
public:
    ElementFem(std::vector<NodeSimple*> nodes, const InterpolationSimple& interpolation)
        : mNodes(nodes)
        , mInterpolation(interpolation)
    {
        assert(mNodes.size() == interpolation.GetNumNodes());
        assert(mNodes.front()->GetNumValues() == interpolation.GetDofDimension());
    }

    ElementFem(std::initializer_list<std::reference_wrapper<NuTo::NodeSimple>> nodes,
               const InterpolationSimple& interpolation)
        : mInterpolation(interpolation)
    {
        for (NuTo::NodeSimple& node : nodes)
            mNodes.push_back(&node);
        assert(mNodes.size() == interpolation.GetNumNodes());
        assert(mNodes.front()->GetNumValues() == interpolation.GetDofDimension());
    }

    virtual NodeValues ExtractNodeValues() const override
    {
        const int dim = GetDofDimension();
        Eigen::VectorXd nodeValues(GetNumNodes() * dim);
        for (size_t i = 0; i < GetNumNodes(); ++i)
            nodeValues.segment(dim * i, dim) = mNodes[i]->GetValues();
        return nodeValues;
    }

    NMatrix GetNMatrix(NaturalCoords ipCoords) const override
    {
        return Matrix::N(Interpolation().GetShapeFunctions(ipCoords), Interpolation().GetNumNodes(),
                         Interpolation().GetDofDimension());
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
        return Interpolation().GetDofDimension();
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

    NodeSimple& GetNode(int i)
    {
        assert(i < mNodes.size());
        return *mNodes[i];
    }


    const NodeSimple& GetNode(int i) const
    {
        assert(i < mNodes.size());
        return *mNodes[i];
    }

private:
    std::vector<NuTo::NodeSimple*> mNodes;
    std::reference_wrapper<const InterpolationSimple> mInterpolation;
};
} /* NuTo */
