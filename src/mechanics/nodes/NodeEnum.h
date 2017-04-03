// $Id$
#pragma once

#include <limits>
#include <map>
#include <set>
#include <string>

namespace NuTo
{
namespace Node
{

enum class eNodeType
{
    NodeCoordinates,
    NodeCoordinatesDof,
    NodeCoordinatesDofNonlocalData,
};


enum class eDof : unsigned char {COORDINATES,
                                 ROTATIONS,
                                 TEMPERATURE,
                                 DISPLACEMENTS,
                                 FINESCALEDISPLACEMENTS,
                                 NONLOCALDATA,
                                 NONLOCALEQPLASTICSTRAIN,
                                 NONLOCALTOTALSTRAIN,
                                 NONLOCALEQSTRAIN,
                                 WATERVOLUMEFRACTION,
                                 RELATIVEHUMIDITY,
                                 ELECTRICPOTENTIAL,
                                 CRACKPHASEFIELD};

//! @brief Gets a set of all Dofs
const std::set<eDof> GetDofSet();




constexpr size_t maxDoFEnumValue = std::numeric_limits<unsigned char>::max();

constexpr unsigned short CombineDofs(eDof rDof1, eDof rDof2)
{
    return static_cast<unsigned short>(rDof1) * maxDoFEnumValue + static_cast<unsigned short>(rDof2);
}

const std::map<eDof, std::string> GetDofMap();


const std::string DofToString(eDof rDof);

eDof DofToEnum(std::string rDof);


}// namespace Node
}// namespace NuTo
