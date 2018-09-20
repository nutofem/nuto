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
    MeshFem(GeometryMeshFem& geometryMesh);

    MeshFem(const MeshFem&) = delete;
    MeshFem& operator=(const MeshFem&) = delete;

    MeshFem(MeshFem&&) = default;
    MeshFem& operator=(MeshFem&&) = default;

    //! @brief adds a clone of `interpolation` to the mesh (prototype pattern)
    //! @param interpolation interpolation that is cloned and added
    //! @return reference to the cloned object
    InterpolationSimple& CreateInterpolation(const InterpolationSimple& interpolation);

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

    //! @brief selects all nodes of `dofType`
    //! @return group containing all selected nodes
    Group<DofNode> NodesTotal(DofType dofType);

    //! @brief selects all element collections
    //! @return group containing all element collections
    Group<ElementCollectionFem> ElementsTotal();

    //! Adds `numInstances` instances of zeros to all nodes of type `dofType`
    //! @param dofType dof type
    void AllocateDofInstances(DofType dofType, int numInstances);

    DofNode& AddNode(Eigen::VectorXd data)
    {
        auto& nd = Nodes.Add(data);
        return nd;
    }

    DofNode& AddNode(double data)
    {
        auto& nd = Nodes.Add(data);
        return nd;
    }

    DofNode& GetNode(int i)
    {
        return Nodes[i];
    }

    size_t NumNodes() const
    {
        return Nodes.Size();
    }

    GeometryMeshFem& GetGeometryMesh() const
    {
        return mGeometryMesh;
    }

    ValueVector<ElementCollectionFem>& GetElements()
    {
        return Elements;
    }

    ElementCollectionFem& GetElement(int i)
    {
        return Elements[i];
    }

    const ElementCollectionFem& GetElement(int i) const
    {
        return Elements[i];
    }

    size_t NumElements() const
    {
        return Elements.Size();
    }

private:
    ValueVector<ElementCollectionFem> Elements;
    GeometryMeshFem& mGeometryMesh;
    ValueVector<DofNode> Nodes;
    std::vector<std::unique_ptr<InterpolationSimple>> mInterpolations;
};
} /* NuTo */
