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
    // Now add the corresponding interpolation
    // rMesh->CreateInterpolation(elm.Interpolation().GetEdgeInterpolation()) -> Has to be implemented for all
    // interpolations
    InterpolationTriangleQuadratic ipol;
    rMesh->CreateInterpolation(ipol);

    // Now add elements
    std::vector<NodeSimple*> nodes0;
    for (int i : ipol.EdgeNodeIds(0))
    {
        nodes0.push_back(&(elm.CoordinateElement().GetNode(i)));
    }
    auto& e0 = rMesh->Elements.Add({{nodes0, *ipol.EdgeInterpolation(0)}});


    std::vector<NodeSimple*> nodes1;
    for (int i : ipol.EdgeNodeIds(1))
    {
        nodes1.push_back(&(elm.CoordinateElement().GetNode(i)));
    }
    auto& e1 = rMesh->Elements.Add({{nodes1, *ipol.EdgeInterpolation(1)}});


    std::vector<NodeSimple*> nodes2;
    for (int i : ipol.EdgeNodeIds(2))
    {
        nodes2.push_back(&(elm.CoordinateElement().GetNode(i)));
    }
    auto& e2 = rMesh->Elements.Add({{nodes2, *ipol.EdgeInterpolation(2)}});

    Group<ElementCollectionFem> edgeElements;
    edgeElements.Add({e0});
    edgeElements.Add({e1});
    edgeElements.Add({e2});
    return edgeElements;
}


} /* NuTo */
