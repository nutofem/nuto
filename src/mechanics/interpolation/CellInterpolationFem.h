#pragma once

#include <vector>
#include "mechanics/nodes/NodeSimple.h"
#include "mechanics/interpolation/CellInterpolationBase.h"
#include "mechanics/interpolation/InterpolationSimple.h"
#include "mechanics/cell/Matrix.h"

namespace NuTo
{

class CellInterpolationFem : public CellInterpolationBase
{
public:
    CellInterpolationFem(std::vector<NuTo::NodeSimple*> nodes, const InterpolationSimple& interpolation)
        : mNodes(nodes)
        , mInterpolation(interpolation)
    {
    }

    virtual NodeValues ExtractNodeValues() const override
    {
        int dim = mNodes[0]->GetNumValues();
        Eigen::VectorXd nodeValues(mNodes.size() * dim);
        for (size_t i = 0; i < mNodes.size(); ++i)
            nodeValues.segment(dim * i, dim) = mNodes[i]->GetValues();
        return nodeValues;
    }

    NMatrix GetNMatrix(NaturalCoords ipCoords) const override
    {
        return Matrix::N(mInterpolation.GetShapeFunctions(ipCoords), mInterpolation.GetNumNodes(),
                         mInterpolation.GetDofDimension());
    }

    Eigen::VectorXd GetShapeFunctions(Eigen::VectorXd ipCoords) const override
    {
        return mInterpolation.GetShapeFunctions(ipCoords);
    }

    Eigen::MatrixXd GetDerivativeShapeFunctions(Eigen::VectorXd ipCoords) const override
    {
        return mInterpolation.GetDerivativeShapeFunctions(ipCoords);
    }

    int GetDofDimension() const override
    {
        return mInterpolation.GetDofDimension();
    }

    int GetNumNodes() const override
    {
        return mInterpolation.GetNumNodes();
    }

private:
    std::vector<NuTo::NodeSimple*> mNodes;
    const InterpolationSimple& mInterpolation;
};
} /* NuTo */
