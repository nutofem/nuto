#pragma once
#include <boost/ptr_container/ptr_vector.hpp>
#include "mechanics/nodes/NodeSimple.h"
#include "mechanics/elements/ElementFem.h"

namespace NuTo
{
class Mesh
{
public:

    NodeSimple& CreateNode(Eigen::VectorXd values)
    {
        mNodes.push_back(new NodeSimple(values));
        return *mNodes.rbegin();
    }

    ElementInterface& CreateElement(std::vector<NodeSimple*> nodes, const InterpolationSimple& interpolation)
    {
        mElements.push_back(new ElementFem{nodes, interpolation});
        return *mElements.rbegin();
    }

    InterpolationSimple& CreateInterpolation(const InterpolationSimple& interpolation)
    {
        mInterpolations.push_back(interpolation.Clone().release());
        return *mInterpolations.rbegin();
    }

private:
    boost::ptr_vector<InterpolationSimple> mInterpolations;
    boost::ptr_vector<ElementFem> mElements;
    boost::ptr_vector<NodeSimple> mNodes;
};
} /* NuTo */
