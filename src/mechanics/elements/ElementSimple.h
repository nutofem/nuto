#pragma once

#include <vector>
#include "mechanics/nodes/NodeSimple.h"
#include "mechanics/interpolation/InterpolationSimple.h"

namespace NuTo
{

class ElementSimple
{
public:
    ElementSimple(std::vector<NuTo::NodeSimple*> rNodes, const InterpolationSimple& rInterpolation)
        : mNodes(rNodes)
        , mInterpolation(rInterpolation)
    {
    }

    //! @brief extracts all node values of this element
    //! @remark virtual to make it testable
    virtual NodeValues ExtractNodeValues() const
    {
        int dim = mNodes[0]->GetNumValues();
        Eigen::VectorXd nodeValues(mNodes.size() * dim);
        for (size_t i = 0; i < mNodes.size(); ++i)
            nodeValues.segment(dim * i, dim) = mNodes[i]->GetValues();
        return nodeValues;
    }

    Eigen::VectorXd Interpolate(const NaturalCoords& rNaturalCoords) const
    {
        return mInterpolation.GetN(rNaturalCoords) * ExtractNodeValues();
    }

    //this has to be removed and replaced by a standard interface to obtain
    //directly shape functions and derivatives
    const InterpolationSimple& GetInterpolation() const
    {
        return mInterpolation;
    }

private:
    std::vector<NuTo::NodeSimple*> mNodes;
    const InterpolationSimple& mInterpolation;
};
} /* NuTo */
