#pragma once
#include <boost/ptr_container/ptr_vector.hpp>
#include "base/ValueVector.h"
#include "mechanics/nodes/NodeSimple.h"
#include "mechanics/elements/ElementCollection.h"


namespace NuTo
{
class MeshFem
{
public:

    InterpolationSimple& CreateInterpolation(const InterpolationSimple& interpolation)
    {
        mInterpolations.push_back(interpolation.Clone().release());
        return *mInterpolations.rbegin();
    }

public:
    ValueVector<NodeSimple> Nodes;
    ValueVector<ElementCollectionFem> Elements;

private:
    boost::ptr_vector<InterpolationSimple> mInterpolations;
};
} /* NuTo */
