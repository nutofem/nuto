#pragma once
#include "nuto/base/Group.h"
#include "nuto/base/ValueVector.h"
#include "nuto/mechanics/DirectionEnum.h"
#include "nuto/mechanics/nodes/DofNode.h"
#include "nuto/mechanics/elements/ElementCollection.h"

#include <memory>
#include <vector>

namespace NuTo
{

class GeometryMeshFem;

//! @brief contains the nodes, elements and interpolations for a classic finite element mesh
//! @remark Elements contain references to nodes. Thus, a copy of MeshFem is not trivially possible and only move is
//! allowed
class MeshFem
{
public:
    MeshFem(const GeometryMeshFem& geometryMesh);

    MeshFem(const MeshFem&) = delete;
    MeshFem& operator=(const MeshFem&) = delete;

    MeshFem(MeshFem&&) = default;
    MeshFem& operator=(MeshFem&&) = default;

    //! @brief adds a clone of `interpolation` to the mesh (prototype pattern)
    //! @param interpolation interpolation that is cloned and added
    //! @return reference to the cloned object
    InterpolationSimple& CreateInterpolation(const InterpolationSimple& interpolation);

    //! @brief selects a coordinate at given `coords`
    //! @param coords global coordinates
    //! @param tol selection tolerance
    //! @return reference to the selected node, throws, if no node is found
    CoordinateNode& NodeAtCoordinate(Eigen::VectorXd coords, double tol = 1.e-10);

    //! @brief selects a node of type `dofType` at given `coords`
    //! @param coords global coordinates
    //! @param dofType dof type
    //! @param tol selection tolerance
    //! @return reference to the selected node, throws, if no node is found
    DofNode& NodeAtCoordinate(Eigen::VectorXd coords, DofType dofType, double tol = 1.e-10);

    //! @brief selects all nodes of type `dofType` where the `coord` in `direction` is within `tol`
    //! @param direction ::X, ::Y, or ::Z
    //! @param dofType dof type
    //! @param axisOffset distance of the node to the axis
    //! @param tol selection tolerance
    //! @return group with selected nodes, the group may be empty if no nodes were found
    Group<DofNode> NodesAtAxis(eDirection direction, DofType dofType, double axisOffset = 0., double tol = 1.e-10);

    //! @brief selects all coordinate nodes where the `coord` in `direction` is within `tol`
    //! @param direction ::X, ::Y, or ::Z
    //! @param axisOffset distance of the node to the axis
    //! @param tol selection tolerance
    //! @return group with selected nodes, the group may be empty if no nodes were found
    Group<CoordinateNode> NodesAtAxis(eDirection direction, double axisOffset = 0., double tol = 1.e-10);

    //! @brief selects all nodes of `dofType`
    //! @return group containing all selected nodes
    Group<DofNode> NodesTotal(DofType dofType);

    //! @brief selects all coordinate nodes
    //! @return group containing all selected nodes
    Group<CoordinateNode> NodesTotal();

    //! @brief selects all element collections
    //! @return group containing all element collections
    Group<ElementCollectionFem> ElementsTotal();

    //! Adds `numInstances` instances of zeros to all nodes of type `dofType`
    //! @param dofType dof type
    void AllocateDofInstances(DofType dofType, int numInstances);

public:
    const GeometryMeshFem& mGeometryMesh;
    ValueVector<DofNode> Nodes;
    ValueVector<ElementCollectionFem> Elements;

private:
    std::vector<std::unique_ptr<InterpolationSimple>> mInterpolations;
};
} /* NuTo */
