#pragma once
#include <boost/ptr_container/ptr_vector.hpp>
#include "base/Group.h"
#include "base/ValueVector.h"
#include "mechanics/DirectionEnum.h"
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

    //! @brief selects a node of type `dofType` at given `coords`
    //! @param coords global coordinates
    //! @param dofType dof type
    //! @param tol selection tolerance
    //! @return reference to the selected node, throws, if no node is found
    NodeSimple& NodeAtCoordinate(Eigen::VectorXd coords, DofType dofType, double tol = 1.e-10)
    {
        for (auto& element : this->Elements)
        {
            auto& dofElement = element.DofElement(dofType);
            const auto& dofInterpolation = dofElement.Interpolation();
            for (int iNode = 0; iNode < dofInterpolation.GetNumNodes(); ++iNode)
            {
                NaturalCoords dofNodeCoords = dofInterpolation.GetLocalCoords(iNode);
                Eigen::VectorXd globalNodeCoords = Interpolate(element.CoordinateElement(), dofNodeCoords);
                if ((globalNodeCoords - coords).isMuchSmallerThan(tol))
                    return dofElement.GetNode(iNode);
            }
        }
        std::stringstream coordsString;
        coordsString << coords.transpose();
        throw NuTo::Exception(__PRETTY_FUNCTION__,
                              "There is no node for dof type " + dofType.GetName() + " at " + coordsString.str());
    }

    //! @brief selects all nodes of type `dofType` where the `coord` in `direction` is within `tol`
    //! @param direction ::X, ::Y, or ::Z
    //! @param dofType dof type
    //! @param axisOffset distance of the node to the axis
    //! @param tol selection tolerance
    //! @return group with selected nodes, the group may be empty if no nodes were found
    Groups::Group<NodeSimple> NodesAtAxis(eDirection direction, DofType dofType, double axisOffset = 0., double tol = 1.e-10)
    {
        Groups::Group<NodeSimple> group;
        const int directionComponent = ToComponentIndex(direction);
        for (auto& element : this->Elements)
        {
            auto& dofElement = element.DofElement(dofType);
            const auto& dofInterpolation = dofElement.Interpolation();
            for (int iNode = 0; iNode < dofInterpolation.GetNumNodes(); ++iNode)
            {
                NaturalCoords dofNodeCoords = dofInterpolation.GetLocalCoords(iNode);
                Eigen::VectorXd globalNodeCoords = Interpolate(element.CoordinateElement(), dofNodeCoords);
                
                if (std::abs(globalNodeCoords[directionComponent] - axisOffset) < tol)
                    group.Add(dofElement.GetNode(iNode));
            }
        }
        return group;
    }

public:
    ValueVector<NodeSimple> Nodes;
    ValueVector<ElementCollectionFem> Elements;

private:
    boost::ptr_vector<InterpolationSimple> mInterpolations;
};


} /* NuTo */
