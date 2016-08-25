#include "ElementDataEnum.h"

#include <boost/algorithm/string.hpp>
#include "nuto/mechanics/MechanicsException.h"


NuTo::ElementData::eElementDataType NuTo::ElementData::ElementDataTypeToEnum(const std::string &rElementDataType)
{
    std::string uppercase = boost::to_upper_copy(rElementDataType);

    ElementData::eElementDataType elementDataType;
    if (uppercase == "CONSTITUTIVELAWIP")
    {
        elementDataType = NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP;
    }
    else if (uppercase == "CONSTITUTIVELAWIPCRACK")
    {
        elementDataType = NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIPCRACK;
    }
    else if (uppercase == "CONSTITUTIVELAWIPNONLOCAL")
    {
        elementDataType = NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIPNONLOCAL;
    }
    else
    {
        throw MechanicsException("[NuTo::ElementData::ElementDataTypeToEnum] Element data type " + rElementDataType + " does not exist.");
    }
    return elementDataType;
}
