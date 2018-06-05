#pragma once
#include "nuto/base/Group.h"
#include "nuto/mechanics/elements/ElementCollection.h"
#include "nuto/mechanics/mesh/MeshFem.h"
#include "nuto/mechanics/interpolation/InterpolationTriangleQuadratic.h"

namespace NuTo
{

//! @brief Adds edge elements to mesh
//! @param rMesh fem mesh, return argument with r and weird pointer syntax to make it clear
//! @param rMesh elm, element whose edges are added
Group<ElementCollectionFem> AddEdgeElements(MeshFem* rMesh, ElementCollectionFem& elm)
{
    auto& elmIpol = elm.CoordinateElement().Interpolation();
    Group<ElementCollectionFem> edgeElements;

    // Now add elements
    for (int edgeId = 0; edgeId < elmIpol.NumEdges(); edgeId++)
    {
        std::vector<NodeSimple*> nodes;
        for (int i : elmIpol.EdgeNodeIds(edgeId))
        {
            nodes.push_back(&(elm.CoordinateElement().GetNode(i)));
        }
        auto& edgeIpol = rMesh->CreateInterpolation(*elmIpol.EdgeInterpolation(edgeId));
        auto& e0 = rMesh->Elements.Add({{nodes, edgeIpol}});
        edgeElements.Add({e0});
    }
    return edgeElements;
}


} /* NuTo */
