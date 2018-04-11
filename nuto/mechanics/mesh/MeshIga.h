#pragma once
#include "nuto/base/Group.h"
#include "nuto/base/ValueVector.h"
#include "nuto/mechanics/DirectionEnum.h"
#include "nuto/mechanics/nodes/NodeSimple.h"
#include "nuto/mechanics/elements/ElementCollection.h"
#include "nuto/mechanics/iga/Nurbs.h"

#include <memory>
#include <vector>

namespace NuTo
{
template <int TDimParameter>
class MeshIga
{
public:
    MeshIga() = default;

    MeshIga(const MeshIga&) = delete;
    MeshIga& operator=(const MeshIga&) = delete;

    MeshIga(MeshIga&&) = default;
    MeshIga& operator=(MeshIga&&) = default;

    Group<NodeSimple> NodesTotal()
    {
        throw Exception(__PRETTY_FUNCTION__, "Iga - Not implemented yet!");
    }

    Group<NodeSimple> NodesTotal(DofType d)
    {
        throw Exception(__PRETTY_FUNCTION__, "Iga - Not implemented yet!");
    }

    Group<ElementCollectionIga<TDimParameter>> ElementsTotal()
    {
        throw Exception(__PRETTY_FUNCTION__, "Iga - Not implemented yet!");
    }

public:
    //! @brief For each dof a nurbs is stored (kind of interpolation)
    DofContainer<Nurbs<TDimParameter>> mDofElements;

    //! @brief Isogeometric elements (spans of knot vector)
    ValueVector<ElementCollectionIga<TDimParameter>> mElements;

    //! @brief Isogeometric elements (spans of knot vector)
    ValueVector<NodeSimple> mNodes;
};
} /* NuTo */
