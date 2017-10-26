#pragma once
#include <boost/ptr_container/ptr_vector.hpp>
#include "base/Group.h"
#include "base/ValueVector.h"
#include "mechanics/DirectionEnum.h"
#include "mechanics/nodes/NodeSimple.h"
#include "mechanics/elements/ElementCollection.h"

namespace NuTo
{
class MeshFem
{
public:
    MeshFem() = default;

    MeshFem(const MeshFem&) = delete;
    MeshFem& operator=(const MeshFem&) = delete;

    MeshFem(MeshFem&&) = default;
    MeshFem& operator=(MeshFem&&) = default;

    InterpolationSimple& CreateInterpolation(const InterpolationSimple& interpolation);


    //! @brief selects a coordinate at given `coords`
    //! @param coords global coordinates
    //! @param tol selection tolerance
    //! @return reference to the selected node, throws, if no node is found
    NodeSimple& NodeAtCoordinate(Eigen::VectorXd coords, double tol = 1.e-10);

    //! @brief selects a node of type `dofType` at given `coords`
    //! @param coords global coordinates
    //! @param dofType dof type
    //! @param tol selection tolerance
    //! @return reference to the selected node, throws, if no node is found
    NodeSimple& NodeAtCoordinate(Eigen::VectorXd coords, DofType dofType, double tol = 1.e-10);

    //! @brief selects all nodes of type `dofType` where the `coord` in `direction` is within `tol`
    //! @param direction ::X, ::Y, or ::Z
    //! @param dofType dof type
    //! @param axisOffset distance of the node to the axis
    //! @param tol selection tolerance
    //! @return group with selected nodes, the group may be empty if no nodes were found
    Groups::Group<NodeSimple> NodesAtAxis(eDirection direction, DofType dofType, double axisOffset = 0.,
                                          double tol = 1.e-10);

    //! @brief selects all coordinate nodes where the `coord` in `direction` is within `tol`
    //! @param direction ::X, ::Y, or ::Z
    //! @param axisOffset distance of the node to the axis
    //! @param tol selection tolerance
    //! @return group with selected nodes, the group may be empty if no nodes were found
    Groups::Group<NodeSimple> NodesAtAxis(eDirection direction, double axisOffset = 0., double tol = 1.e-10);

    //! @brief selects all nodes of `dofType`
    //! @return group containing all selected nodes
    Groups::Group<NodeSimple> NodesTotal(DofType dofType);

    //! @brief selects all coordinate nodes
    //! @return group containing all selected nodes
    Groups::Group<NodeSimple> NodesTotal();


public:
    ValueVector<NodeSimple> Nodes;
    ValueVector<ElementCollectionFem> Elements;

private:
    boost::ptr_vector<InterpolationSimple> mInterpolations;
};
} /* NuTo */
