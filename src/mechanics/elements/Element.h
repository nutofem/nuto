#pragma once

#include <vector>
#include "mechanics/nodes/NodeSimple.h"
#include "mechanics/interpolation/Interpolation.h"

namespace NuTo
{

class Element
{
public:
    Element(std::vector<NuTo::NodeSimple*> rNodes, const Interpolation& rInterpolation)
        : mNodes(rNodes)
        , mInterpolation(rInterpolation)
    {
    }

    Eigen::VectorXd ExtractNodeValues() const
    {
        int dim = mNodes[0]->GetNumValues();
        Eigen::VectorXd nodeValues(mNodes.size() * dim);
        for (size_t i = 0; i < mNodes.size(); ++i)
            nodeValues.segment(dim * i, dim) = mNodes[i]->GetValues();
        return nodeValues;
    }

    Eigen::VectorXd Interpolate(const Eigen::VectorXd& rLocalCoords) const
    {
        return GetN(rLocalCoords) * ExtractNodeValues();
    }

    const Interpolation& GetInterpolation() const
    {
        return mInterpolation;
    }

private:
    //! @brief 'blows' shape functions to N matrix
    Eigen::MatrixXd GetN(const Eigen::VectorXd& rLocalCoords) const
    {
        int dim = mNodes[0]->GetNumValues();
        Eigen::MatrixXd N(dim, dim * mNodes.size());

        auto shapeFunctions = mInterpolation.GetShapeFunctions(rLocalCoords);

        for (size_t i = 0; i < mNodes.size(); ++i)
            N.block(0, i * dim, dim, dim) = Eigen::MatrixXd::Identity(dim, dim) * shapeFunctions[i];
        return N;
    }


    std::vector<NuTo::NodeSimple*> mNodes;
    const Interpolation& mInterpolation;
};
} /* NuTo */
