#pragma once
#include <vector>
#include <boost/ptr_container/ptr_vector.hpp>
#include "mechanics/elements/ElementSimple.h"
#include "mechanics/nodes/NodeSimple.h"
#include "mechanics/interpolation/InterpolationSimple.h"

namespace NuTo
{
class Mesh
{
public:
    Mesh()
    {
        mNodes.reserve(1e6);
        mElements.reserve(1e6);
    }

    NodeSimple& CreateNode(Eigen::VectorXd rValues)
    {
        mNodes.push_back(NodeSimple(rValues));
        return *mNodes.rbegin();
    }

    ElementSimple& CreateElement(std::vector<NodeSimple*> rNodes, const InterpolationSimple& rInterpolation)
    {
        mElements.push_back({rNodes, rInterpolation});
        return *mElements.rbegin(); 
    }

    InterpolationSimple& CreateInterpolation(const InterpolationSimple& rInterpolation)
    {
        mInterpolations.push_back(rInterpolation.Clone().release());
        return *mInterpolations.rbegin();
    }

private:
    boost::ptr_vector<InterpolationSimple> mInterpolations;
    std::vector<ElementSimple> mElements;
    std::vector<NodeSimple> mNodes;
};
} /* NuTo */
