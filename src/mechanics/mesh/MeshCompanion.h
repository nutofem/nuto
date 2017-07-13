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
//! @param rS almighty Structure
//! @param rFileName File name
//! @return Vector of pair, with element.first containing the group id, and element.second the interpolation type id

std::vector<std::pair<int, int>> ImportFromGmsh(Structure& rS, const std::string &rFileName);
std::vector<std::pair<int, int>> ImportFromGmsh(Structure& rS, const std::string &rFileName, bool useNewNumbers);
std::vector<std::pair<int, int>> ImportFromGmsh(Structure& rS, const std::string &rFileName, bool useNewNumbers, std::map<int, int> &newNodes, std::map<int, int> &gmshNodes);


//! @brief changes the node structure to match the interpolation type for all elements
//! the node merge distance and the box size are calculated from the element sizes
//! @param rS almighty Structure
void ElementTotalConvertToInterpolationType(Structure& rS);

//! @brief changes the node structure to match the interpolation type
//! the node merge distance and the box size are calculated from the element sizes
//! @param rS almighty Structure
//! @param rGroupNumberElements group for elements (coordinates only) to be converted
void ElementConvertToInterpolationType(Structure& rS, int rGroupNumberElements);

//! @brief changes the node structure to match the interpolation type for all elements
//! @param rS almighty Structure
//! @param rNodeDistanceMerge Distance of nodes to be joined (should be significantly smaller than the node distance in
//! the mesh)
void ElementTotalConvertToInterpolationType(Structure& rS, double rNodeDistanceMerge);

//! @brief changes the node structure to match the interpolation type
//! @param rS almighty Structure
//! @param rGroupNumberElements group for elements (coordinates only) to be converted
//! @param rNodeDistanceMerge Distance of nodes to be joined (should be significantly smaller than the node distance in
//! the mesh)
void ElementConvertToInterpolationType(Structure& rS, int rGroupNumberElements, double rNodeDistanceMerge);

//! @brief creates PRISM3D elements between two element groups
//! @param rS almighty Structure (3D!)
//! @param rGroupMaster group of tetrahedron elements on one side of the prisms
//! @param rGroupSlave group of tetrahedron elements on the other side of the prisms
//! @param rThickness thickness of the newly created PRISM3D elements (height of the prisms)
//! @return pair<PrismGroupId, PrismInterpolationTypeId> (like in ImportFromGmsh)
std::pair<int, int> ElementPrismsCreate(NuTo::Structure& rS, int rGroupMaster, int rGroupSlave, double rThickness);


//! @brief builds a vector of all element pointers inside an element group
//! @param rS almighty Structure
//! @param rGroupNumberElements element group id
//! @return vector of all element pointers inside an element group
std::vector<ElementBase*> GetElementVector(Structure& rS, int rElementGroupId);

//! @brief builds a mapping from node ptr to node id
//! @param rS almighty Structure
//! @return mapping from node ptr to node id
std::map<const NuTo::NodeBase*, int> GetNodeToIdMap(Structure& rS);

} // namespace MeshCompanion
} // namespace NuTo
