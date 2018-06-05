#pragma once
#include "nuto/base/Group.h"
#include "nuto/mechanics/elements/ElementCollection.h"
#include "nuto/mechanics/mesh/MeshFem.h"
#include "nuto/mechanics/interpolation/InterpolationTriangleQuadratic.h"
#include <set>

namespace NuTo
{

//! @brief Adds edge elements to mesh
//! @param rMesh fem mesh, return argument with r and weird pointer syntax to make it clear
//! @param rMesh elm, element whose edges are added
Group<ElementCollectionFem> AddEdgeElements(MeshFem* rMesh, ElementCollectionFem& elmColl)
{
    Group<ElementCollectionFem> edgeElements;

    auto& elmIpol = elmColl.CoordinateElement().Interpolation();
    for (int edgeId = 0; edgeId < elmIpol.NumEdges(); edgeId++)
    {
        std::vector<NodeSimple*> nodes;
        for (int i : elmIpol.EdgeNodeIds(edgeId))
        {
            nodes.push_back(&(elmColl.CoordinateElement().GetNode(i)));
        }
        auto& edgeIpol = rMesh->CreateInterpolation(*elmIpol.EdgeInterpolation(edgeId));
        auto& e0 = rMesh->Elements.Add({{nodes, edgeIpol}});
        edgeElements.Add({e0});
    }
    return edgeElements;
}

//! @brief Adds edge elements to mesh
//! @param rMesh fem mesh, return argument with r and weird pointer syntax to make it clear
//! @param rMesh elms, elements whose edges are added, no duplicates
//! edge direction is from smallest pointer to largest
Group<ElementCollectionFem> AddEdgeElements(MeshFem* rMesh, Group<ElementCollectionFem>& elmCollGroup)
{
    Group<ElementCollectionFem> edgeElements;

    auto cmp = [](std::vector<NodeSimple*> a, std::vector<NodeSimple*> b) {
        for (int i = 0; i < a.size(); i++)
        {
            if (a[i] < b[i])
                return true;
            else if (a[i] > b[i])
            {
                return false;
            }
        }
        return false;
    };
    std::set<std::vector<NodeSimple*>, decltype(cmp)> alreadyKnownEdges(cmp);

    for (ElementCollectionFem& elmColl : elmCollGroup)
    {
        auto& elmIpol = elmColl.CoordinateElement().Interpolation();
        for (int edgeId = 0; edgeId < elmIpol.NumEdges(); edgeId++)
        {
            std::vector<NodeSimple*> nodes;
            for (int i : elmIpol.EdgeNodeIds(edgeId))
            {
                nodes.push_back(&(elmColl.CoordinateElement().GetNode(i)));
            }
            if (nodes.front() < nodes.back())
            {
                std::reverse(nodes.begin(), nodes.end());
            }
            if (alreadyKnownEdges.insert(nodes).second)
            {
                auto& edgeIpol = rMesh->CreateInterpolation(*elmIpol.EdgeInterpolation(edgeId));
                auto& e0 = rMesh->Elements.Add({{nodes, edgeIpol}});
                edgeElements.Add({e0});
            }
        }
    }
    return edgeElements;
}


} /* NuTo */
