// $Id$
#ifndef ELEMENTDATAENUM_H_
#define ELEMENTDATAENUM_H_

#include <boost/algorithm/string.hpp>
#include "nuto/mechanics/MechanicsException.h"

namespace NuTo
{
namespace ElementData
{
enum eElementDataType
{
    NOELEMENTDATA = 0, CONSTITUTIVELAWIP,                      //!< constitutive law and integration points
    CONSTITUTIVELAWIPNONLOCAL,              //!< constitutive law, integration points and nonlocal data
    CONSTITUTIVELAWIPCRACK,                  //!< constitutive law, integration points and crack data
    VARIABLECONSTITUTIVELAWIP,              //!< variable constitutive laws and integration points
};

static inline eElementDataType ElementDataTypeToEnum(const std::string& rElementDataType)
{
    std::string uppercase = boost::to_upper_copy(rElementDataType);

    ElementData::eElementDataType elementDataType;
    if (uppercase == "CONSTITUTIVELAWIP")
    {
        elementDataType = NuTo::ElementData::CONSTITUTIVELAWIP;
    }
    else if (uppercase == "CONSTITUTIVELAWIPCRACK")
    {
        elementDataType = NuTo::ElementData::CONSTITUTIVELAWIPCRACK;
    }
    else if (uppercase == "CONSTITUTIVELAWIPNONLOCAL")
    {
        elementDataType = NuTo::ElementData::CONSTITUTIVELAWIPNONLOCAL;
    }
    else
    {
        throw MechanicsException("[NuTo::ElementData::ElementDataTypeToEnum] Element data type " + rElementDataType + " does not exist.");
    }
    return elementDataType;
}
}
}
#endif /* ELEMENTDATAENUM_H_ */
