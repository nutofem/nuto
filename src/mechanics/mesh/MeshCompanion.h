#pragma once

#include <vector>
#include <string>


namespace NuTo
{

class Structure;

namespace MeshCompanion
{

//! @brief Import from gmsh
//!        Creates groups according to gmsh's physical entities and creates an interpolation types for each group
//! @param rS almighty Structure
//! @param rFileName File name
//! @return Vector of pair, with element.first containing the group id, and element.second the interpolation type id
std::vector<std::pair<int, int>> ImportFromGmsh(Structure &rS, const std::string &rFileName);

} // namespace MeshCompanion
} // namespace NuTo