#pragma once
#include <boost/ptr_container/ptr_vector.hpp>
#include "mechanics/interpolation/CellInterpolationBase.h"
#include "mechanics/nodes/NodeSimple.h"
#include "mechanics/interpolation/CellInterpolationFem.h"

namespace NuTo
{
class Mesh
{
public:

    NodeSimple& CreateNode(Eigen::VectorXd rValues)
    {
        mNodes.push_back(new NodeSimple(rValues));
        return *mNodes.rbegin();
    }

    CellInterpolationBase& CreateElement(std::vector<NodeSimple*> rNodes, const InterpolationSimple& rInterpolation)
    {
        mElements.push_back(new CellInterpolationFem{rNodes, rInterpolation});
        return *mElements.rbegin();
    }

    InterpolationSimple& CreateInterpolation(const InterpolationSimple& rInterpolation)
    {
        mInterpolations.push_back(rInterpolation.Clone().release());
        return *mInterpolations.rbegin();
    }

private:
    boost::ptr_vector<InterpolationSimple> mInterpolations;
    boost::ptr_vector<CellInterpolationFem> mElements;
    boost::ptr_vector<NodeSimple> mNodes;
};
} /* NuTo */
