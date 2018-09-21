#pragma once
#include "nuto/base/Group.h"
#include "nuto/mechanics/elements/ElementCollection.h"
#include "nuto/mechanics/mesh/MeshFem.h"

namespace NuTo
{

//! @brief Adds edge elements to mesh if addEdge is true
//! @param rMesh fem mesh, return argument with r and weird pointer syntax to make it clear
//! @param elmCollGroup, elements whose edges are added
//! @param addEdge, function that takes an edge and returns true if it should be included
//! @param orientedEdges if set to false, edges with reversed direction are considered equal
//! (edge nodes are sorted by pointer so that the first is smaller than the last)
Group<ElementCollectionFem> AddEdgeElements(MeshFem* rMesh, Group<ElementCollectionFem>& elmCollGroup,
                                            std::function<bool(ElementCollectionFem&)> addEdge,
                                            bool orientedEdges = false);

//! @brief Adds edge elements to mesh (coordinate elements)
//! @param rMesh fem mesh, return argument with r and weird pointer syntax to make it clear
//! @param elmCollGroup, elements whose edges are added
//! @param addFace, function that takes a face and returns true if it should be included
//! @param orientedEdges if set to false, edges with reversed direction are considered equal
//! (edge nodes are sorted by pointer so that the first is smaller than the last)
Group<ElementCollectionFem> AddEdgeElements(MeshFem* rMesh, Group<ElementCollectionFem>& elmCollGroup,
                                            bool orientedEdges = false);

//! @brief Adds face elements to mesh if addFace is true
//! @param rMesh fem mesh, return argument with r and weird pointer syntax to make it clear
//! @param elmCollGroup, elements whose faces are added
//! @paragraph orientedFaces if set to false, faces with reversed orientation are considered equal
//! (face nodes are sorted by pointer so that the first is smallest. If orientedFaces is set to false
//! then the orientation is chosen that gives the smallest second node)
Group<ElementCollectionFem> AddFaceElements(MeshFem* rMesh, Group<ElementCollectionFem>& elmCollGroup,
                                            std::function<bool(ElementCollectionFem&)> addFace,
                                            bool orientedFaces = false);

//! @brief Adds face elements to mesh
//! @param rMesh fem mesh, return argument with r and weird pointer syntax to make it clear
//! @param elmCollGroup, elements whose faces are added
//! @paragraph orientedFaces if set to false, faces with reversed orientation are considered equal
//! (face nodes are sorted by pointer so that the first is smallest. If orientedFaces is set to false
//! then the orientation is chosen that gives the smallest second node)
Group<ElementCollectionFem> AddFaceElements(MeshFem* rMesh, Group<ElementCollectionFem>& elmCollGroup,
                                            bool orientedFaces = false);

} /* NuTo */
