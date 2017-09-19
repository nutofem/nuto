#pragma once

#include <vector>
#include "mechanics/nodes/NodeSimple.h"
#include "InterpolationSimple.h"
#include "CellInterpolationBase.h"

namespace NuTo
{

class CellInterpolationFEM : public CellInterpolationBase
{
public:
    CellInterpolationFEM(std::vector<NuTo::NodeSimple*> rNodes, const InterpolationSimple& rInterpolation)
        : mNodes(rNodes)
        , mInterpolation(rInterpolation)
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
