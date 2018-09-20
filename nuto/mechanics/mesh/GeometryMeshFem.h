#pragma once
#include "nuto/base/Group.h"
#include "nuto/base/ValueVector.h"
#include "nuto/mechanics/DirectionEnum.h"
#include "nuto/mechanics/nodes/DofNode.h"
#include "nuto/mechanics/elements/ElementCollection.h"
#include "nuto/mechanics/elements/CoordinateElementFem.h"

#include <memory>
#include <vector>

namespace NuTo
{
//! @brief contains the nodes, elements and interpolations for a classic finite element mesh
//! @remark Elements contain references to nodes. Thus, a copy of GeometryMeshFem is not trivially possible and only
//! move is
//! allowed
class GeometryMeshFem
{
public:
    GeometryMeshFem() = default;

    GeometryMeshFem(const GeometryMeshFem&) = delete;
    GeometryMeshFem& operator=(const GeometryMeshFem&) = delete;

    GeometryMeshFem(GeometryMeshFem&&) = default;
    GeometryMeshFem& operator=(GeometryMeshFem&&) = default;

    //! @brief adds a clone of `interpolation` to the mesh (prototype pattern)
    //! @param interpolation interpolation that is cloned and added
    //! @return reference to the cloned object
    InterpolationSimple& CreateInterpolation(const InterpolationSimple& interpolation);

    //! @brief selects a coordinate at given `coords`
    //! @param coords global coordinates
    //! @param tol selection tolerance
    //! @return reference to the selected node, throws, if no node is found
    CoordinateNode& NodeAtCoordinate(Eigen::VectorXd coords, double tol = 1.e-10);

    //! @brief selects all coordinate nodes where the `coord` in `direction` is within `tol`
    //! @param direction ::X, ::Y, or ::Z
    //! @param axisOffset distance of the node to the axis
    //! @param tol selection tolerance
    //! @return group with selected nodes, the group may be empty if no nodes were found
    Group<CoordinateNode> NodesAtAxis(eDirection direction, double axisOffset = 0., double tol = 1.e-10);

    //! @brief selects all coordinate nodes
    //! @return group containing all selected nodes
    Group<CoordinateNode> NodesTotal();

    //! @brief selects all element collections
    //! @return group containing all element collections
    Group<const CoordinateElementFem> ElementsTotal();

    CoordinateNode& AddNode(Eigen::VectorXd data)
    {
        auto& nd = CoordinateNodes.Add(data);
        return nd;
    }

    CoordinateNode& AddNode(double data)
    {
        auto& nd = CoordinateNodes.Add(data);
        return nd;
    }

    CoordinateNode& GetNode(int i)
    {
        return CoordinateNodes[i];
    }

    size_t NumNodes() const
    {
        return CoordinateNodes.Size();
    }

    const CoordinateElementFem& GetElement(int i)
    {
        return Elements[i];
    }

    const CoordinateElementFem& AddElement(std::vector<CoordinateNode*> nodes, const InterpolationSimple& interpolation)
    {
        auto& elm = Elements.Add({nodes, interpolation});
        elm.SetId(Elements.Size() - 1);
        return elm;
    }

    const CoordinateElementFem& AddElement(std::initializer_list<std::reference_wrapper<CoordinateNode>> nodes,
                                           const InterpolationSimple& interpolation)
    {
        auto& elm = Elements.Add({nodes, interpolation});
        elm.SetId(Elements.Size() - 1);
        return elm;
    }

    size_t NumElements() const
    {
        return Elements.Size();
    }

private:
    ValueVector<CoordinateNode> CoordinateNodes;
    ValueVector<CoordinateElementFem> Elements;
    std::vector<std::unique_ptr<InterpolationSimple>> mInterpolations;
};
} /* NuTo */
