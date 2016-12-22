#include "mechanics/elements/ElementEnum.h"


#include <algorithm>
#include "mechanics/MechanicsException.h"

const std::map<NuTo::Element::eElementType, std::string> NuTo::Element::GetElementTypeMap()
{
    const std::map<eElementType, std::string> map =
       {{eElementType::CONTINUUMELEMENT,"CONTINUUMELEMENT"},
        {eElementType::CONTINUUMELEMENTIGA,"CONTINUUMELEMENTIGA"},
        {eElementType::CONTINUUMBOUNDARYELEMENT, "CONTINUUMBOUNDARYELEMENT"},
        {eElementType::CONTINUUMBOUNDARYELEMENTCONSTRAINEDCONTROLNODE, "CONTINUUMBOUNDARYELEMENTCONSTRAINEDCONTROLNODE"},
        {eElementType::ELEMENT1DINXD, "ELEMENT1DINXD"},
        {eElementType::ELEMENT2DINTERFACE, "ELEMENTINTERFACE"}};
    return map;
}

std::string NuTo::Element::ElementTypeToString(NuTo::Element::eElementType rType)
{
    try
    {
        return GetElementTypeMap().find(rType)->second;
    }
    catch (const std::out_of_range& e)
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Enum undefined or not implemented.");
    }
}

NuTo::Element::eElementType NuTo::Element::ElementTypeToEnum(std::string rType)
{
    std::transform(rType.begin(), rType.end(),rType.begin(), ::toupper);

    for(auto entry : GetElementTypeMap())
        if (entry.second == rType)
            return entry.first;

    throw MechanicsException(__PRETTY_FUNCTION__, "ElementType " + rType + " has no enum equivalent or is not implemented.");
}
