#pragma once
#include "nuto/base/Group.h"
#include "nuto/mechanics/elements/ElementCollection.h"
#include "nuto/mechanics/mesh/MeshFem.h"
#include <set>

namespace NuTo
{

//! @brief Adds edge elements to mesh (coordinate elements)
//! @param rMesh fem mesh, return argument with r and weird pointer syntax to make it clear
//! @param elmCollGroup, elements whose edges are added
//! @param orientedEdges if set to false, edges with reversed direction are considered equal
//! (edge nodes are sorted by pointer so that the first is smaller than the last)
Group<ElementCollectionFem> AddEdgeElements(MeshFem* rMesh, Group<ElementCollectionFem>& elmCollGroup,
                                            bool orientedEdges = false)
{
    Group<ElementCollectionFem> edgeElements;

    auto cmp = [](std::vector<NodeSimple*> a, std::vector<NodeSimple*> b) {
        for (size_t i = 0; i < a.size(); i++)
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
            if (!orientedEdges)
            {
                if (nodes.front() < nodes.back())
                {
                    std::reverse(nodes.begin(), nodes.end());
                }
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

//! @brief Adds face elements to mesh
//! @param rMesh fem mesh, return argument with r and weird pointer syntax to make it clear
//! @param elmCollGroup, elements whose faces are added
//! @paragraph orientedFaces if set to false, faces with reversed orientation are considered equal
//! (face nodes are sorted by pointer so that the first is smallest. orientedFaces if set to false
//! then the orientation is chosen that gives the smallest second node)
Group<ElementCollectionFem> AddFaceElements(MeshFem* rMesh, Group<ElementCollectionFem>& elmCollGroup,
                                            bool orientedFaces = false)
{
    Group<ElementCollectionFem> faceElements;

    auto cmp = [](std::vector<NodeSimple*> a, std::vector<NodeSimple*> b) {
        for (size_t i = 0; i < a.size(); i++)
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
    std::set<std::vector<NodeSimple*>, decltype(cmp)> alreadyKnownFaces(cmp);

    for (ElementCollectionFem& elmColl : elmCollGroup)
    {
        auto& elmIpol = elmColl.CoordinateElement().Interpolation();
        for (int faceId = 0; faceId < elmIpol.NumFaces(); faceId++)
        {
            std::vector<NodeSimple*> nodes;
            for (int i : elmIpol.FaceNodeIds(faceId))
            {
                nodes.push_back(&(elmColl.CoordinateElement().GetNode(i)));
            }
            // Sorting nodes circular smallest first
            auto minElm = std::min_element(nodes.begin(), nodes.end());
            std::rotate(nodes.begin(), minElm, nodes.end());
            if (!orientedFaces)
            {
                // chose orientation so that the next smallest is second
                if (nodes[1] > nodes[nodes.size() - 1])
                {
                    std::reverse(nodes.begin(), nodes.end());
                    std::rotate(nodes.begin(), std::min_element(nodes.begin(), nodes.end()), nodes.end());
                }
            }
            if (alreadyKnownFaces.insert(nodes).second)
            {
                auto& faceIpol = rMesh->CreateInterpolation(*elmIpol.FaceInterpolation(faceId));
                auto& e0 = rMesh->Elements.Add({{nodes, faceIpol}});
                faceElements.Add({e0});
            }
        }
    }
    return faceElements;
}

} /* NuTo */
