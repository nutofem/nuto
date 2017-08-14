#pragma once

#include <limits>
#include <set>
#include <string>

namespace NuTo
{
namespace Node
{

enum class eDof : unsigned char
{
    COORDINATES,
    TEMPERATURE,
    DISPLACEMENTS,
    NONLOCALEQSTRAIN,
    WATERVOLUMEFRACTION,
    RELATIVEHUMIDITY,
    ELECTRICPOTENTIAL,
    CRACKPHASEFIELD
};

//! @brief Gets a set of all Dofs
std::set<eDof> GetDofSet();

//! @brief Get number of components of a DOF type
//! @param dofType DOF type
//! @param dimension The dimension of the problem
int GetNumComponents(eDof dofType, int dimension);

constexpr unsigned short CombineDofs(eDof rDof1, eDof rDof2)
{
    return static_cast<unsigned short>(rDof1) * std::numeric_limits<unsigned char>::max() +
           static_cast<unsigned short>(rDof2);
}

std::string DofToString(eDof rDof);

eDof DofToEnum(std::string rDof);

} // namespace Node
} // namespace NuTo
