#pragma once
#include <boost/ptr_container/ptr_vector.hpp>
#include "base/ValueVector.h"
#include "mechanics/nodes/NodeSimple.h"
#include "mechanics/elements/ElementCollection.h"

#include <sstream>

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

const NodeSimple& GetNodeAt(const ValueVector<ElementCollectionFem>& elements, Eigen::VectorXd coords, DofType dofType, double tol = 1.e-10)
{
    for (const auto& element : elements)
    {
        const auto& dofElement = element.DofElement(dofType);
        const auto& dofInterpolation = dofElement.Interpolation();
        for (int iNode = 0; iNode < dofInterpolation.GetNumNodes(); ++iNode)
        {
            NaturalCoords dofNodeCoords = dofInterpolation.GetLocalCoords(iNode);
            Eigen::VectorXd globalNodeCoords = Interpolate(element.CoordinateElement(), dofNodeCoords); 
            if ( (globalNodeCoords - coords).isMuchSmallerThan(tol) )
                return dofElement.GetNode(iNode);
        }
    }
    std::stringstream coordsString;
    coordsString << coords.transpose();
    throw NuTo::Exception(__PRETTY_FUNCTION__, "There is no node for dof type " + dofType.GetName() + " at " + coordsString.str()); 
}

} /* NuTo */
