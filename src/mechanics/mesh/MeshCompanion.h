#pragma once

#include <vector>
#include <string>
#include <map>


namespace NuTo
{

class Structure;
class ElementBase;
class NodeBase;

namespace MeshCompanion
{


//! @brief Import from gmsh
//!        Creates groups according to gmsh's physical entities and creates an interpolation types for each group
//! @param s almighty Structure
//! @param fileName File name
//! @return Vector of pair, with element.first containing the group id, and element.second the interpolation type id
std::vector<std::pair<int, int>> ImportFromGmsh(Structure& s, const std::string& fileName);

//! @brief changes the node structure to match the interpolation type for all elements
//! the node merge distance and the box size are calculated from the element sizes
//! @param s almighty Structure
void ElementTotalConvertToInterpolationType(Structure& s);

//! @brief changes the node structure to match the interpolation type
//! the node merge distance and the box size are calculated from the element sizes
//! @param s almighty Structure
//! @param elementGroupId group for elements (coordinates only) to be converted
void ElementConvertToInterpolationType(Structure& s, int elementGroupId);

//! @brief changes the node structure to match the interpolation type for all elements
//! @param s almighty Structure
//! @param nodeMergeDistance Distance of nodes to be joined (should be significantly smaller than the node distance in
//! the mesh)
void ElementTotalConvertToInterpolationType(Structure& s, double nodeMergeDistance);

//! @brief changes the node structure to match the interpolation type
//! @param s almighty Structure
//! @param elementGroupId group for elements (coordinates only) to be converted
//! @param nodeMergeDistance Distance of nodes to be joined (should be significantly smaller than the node distance in
//! the mesh)
void ElementConvertToInterpolationType(Structure& s, int elementGroupId, double nodeMergeDistance);

//! @brief creates PRISM3D elements between two element groups
//! @param s almighty Structure (3D!)
//! @param groupIdMaster group of tetrahedron elements on one side of the prisms
//! @param groupIdSlave group of tetrahedron elements on the other side of the prisms
//! @param thickness thickness of the newly created PRISM3D elements (height of the prisms)
//! @return pair<PrismGroupId, PrismInterpolationTypeId> (like in ImportFromGmsh)
std::pair<int, int> ElementPrismsCreate(NuTo::Structure& s, int groupIdMaster, int groupIdSlave, double thickness);


//! @brief builds a vector of all element pointers inside an element group
//! @param s almighty Structure
//! @param elementGroupId element group id
//! @return vector of all element pointers inside an element group
std::vector<ElementBase*> GetElementVector(Structure& s, int elementGroupId);

//! @brief builds a mapping from node ptr to node id
//! @param s almighty Structure
//! @return mapping from node ptr to node id
std::map<const NuTo::NodeBase*, int> GetNodeToIdMap(Structure& s);

} // namespace MeshCompanion
} // namespace NuTo
