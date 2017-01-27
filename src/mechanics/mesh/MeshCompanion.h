#pragma once

#include <vector>
#include <string>


namespace NuTo
{

class Structure;
class ElementBase;

namespace MeshCompanion
{

//! @brief Import from gmsh
//!        Creates groups according to gmsh's physical entities and creates an interpolation types for each group
//! @param rS almighty Structure
//! @param rFileName File name
//! @return Vector of pair, with element.first containing the group id, and element.second the interpolation type id
std::vector<std::pair<int, int>> ImportFromGmsh(Structure& rS, const std::string &rFileName);

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
//! @param rNodeDistanceMerge Distance of nodes to be joined (should be significantly smaller than the node distance in the mesh)
void ElementTotalConvertToInterpolationType(Structure& rS, double rNodeDistanceMerge);

//! @brief changes the node structure to match the interpolation type
//! @param rS almighty Structure
//! @param rGroupNumberElements group for elements (coordinates only) to be converted
//! @param rNodeDistanceMerge Distance of nodes to be joined (should be significantly smaller than the node distance in the mesh)
void ElementConvertToInterpolationType(Structure& rS, int rGroupNumberElements, double rNodeDistanceMerge);

std::vector<ElementBase*> GetElementVector(Structure& rS, int rElementGroupId);

} // namespace MeshCompanion
} // namespace NuTo