#pragma once

#include <string>

namespace NuTo
{
namespace ElementData
{
enum class eElementDataType
{
    NOELEMENTDATA = 0,
    CONSTITUTIVELAWIP,                      //!< constitutive law and integration points
    CONSTITUTIVELAWIPNONLOCAL,              //!< constitutive law, integration points and nonlocal data
    VARIABLECONSTITUTIVELAWIP,              //!< variable constitutive laws and integration points
};

eElementDataType ElementDataTypeToEnum(const std::string& rElementDataType);

}// namespace ElementData
}// namespace NuTo
