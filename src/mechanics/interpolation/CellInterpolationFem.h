#pragma once

#include <vector>
#include "mechanics/nodes/NodeSimple.h"
#include "mechanics/interpolation/CellInterpolationBase.h"
#include "mechanics/interpolation/InterpolationSimple.h"
#include "mechanics/cell/Matrix.h"
#include "math/NaturalCoordinateMemoizer.h"

namespace NuTo
{

class CellInterpolationFem : public CellInterpolationBase
{
public:
    CellInterpolationFem(std::vector<NuTo::NodeSimple*> rNodes, const InterpolationSimple& rInterpolation)
        : mNodes(rNodes)
        , mInterpolation(rInterpolation)
        , mShapeFunctions([=](const Eigen::VectorXd& v) { return mInterpolation.GetShapeFunctions(v); })
        , mMatrixN([=](const Eigen::VectorXd& v) {
            return Matrix::N(mInterpolation.GetShapeFunctions(v), GetNumNodes(), GetDofDimension());
        })
        , mDerivativeShapeFunctionsNatural(
                  [=](const Eigen::VectorXd& v) { return mInterpolation.GetDerivativeShapeFunctions(v); })
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
        return mMatrixN.Get(ipCoords);
    }

    Eigen::VectorXd GetShapeFunctions(Eigen::VectorXd ipCoords) const override
    {
        return mShapeFunctions.Get(ipCoords); 
    }

    Eigen::MatrixXd GetDerivativeShapeFunctions(Eigen::VectorXd ipCoords) const override
    {
        return mDerivativeShapeFunctionsNatural.Get(ipCoords); 
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

    template <typename TResult>
    using Memoizer = NuTo::NaturalCoordinateMemoizerMap<TResult, Eigen::VectorXd>;

    Memoizer<Eigen::VectorXd> mShapeFunctions;
    Memoizer<Eigen::MatrixXd> mMatrixN;
    Memoizer<Eigen::MatrixXd> mDerivativeShapeFunctionsNatural;
};
} /* NuTo */
