#include "IpDataEnum.h"

#include <boost/algorithm/string.hpp>
#include "base/Exception.h"


const std::map<NuTo::IpData::eIpStaticDataType, std::string> NuTo::IpData::GetIpStaticDataTypeMap()
{
    const std::map<eIpStaticDataType, std::string> shapeTypeMap = {
            {eIpStaticDataType::DAMAGE, "DAMAGE"},
            {eIpStaticDataType::ELASTIC_ENERGY, "ELASTIC_ENERGY"},
            {eIpStaticDataType::ENGINEERING_PLASTIC_STRAIN, "ENGINEERING_PLASTIC_STRAIN"},
            {eIpStaticDataType::ENGINEERING_STRAIN, "ENGINEERING_STRAIN"},
            {eIpStaticDataType::ENGINEERING_STRESS, "ENGINEERING_STRESS"},
            {eIpStaticDataType::EXTRAPOLATION_ERROR, "EXTRAPOLATION_ERROR"},
            {eIpStaticDataType::HEAT_FLUX, "HEAT_FLUX"},
            {eIpStaticDataType::ELECTRIC_DISPLACEMENT, "ELECTRIC_DISPLACEMENT"},
            {eIpStaticDataType::ELECTRIC_FIELD, "ELECTRIC_FIELD"},
            {eIpStaticDataType::INTERNAL_ENERGY, "INTERNAL_ENERGY"},
            {eIpStaticDataType::LATTICE_PLASTIC_STRAIN, "LATTICE_PLASTIC_STRAIN"},
            {eIpStaticDataType::LATTICE_STRAIN, "LATTICE_STRAIN"},
            {eIpStaticDataType::LATTICE_STRESS, "LATTICE_STRESS"},
            {eIpStaticDataType::SHRINKAGE_STRAIN, "SHRINKAGE_STRAIN"},
            {eIpStaticDataType::THERMAL_STRAIN, "THERMAL_STRAIN"},
            {eIpStaticDataType::TOTAL_INELASTIC_EQ_STRAIN, "TOTAL_INELASTIC_EQUIVALENT_STRAIN"}};
    return shapeTypeMap;
}


std::string NuTo::IpData::IpStaticDataTypeToString(const NuTo::IpData::eIpStaticDataType& rIpStaticDataType)
{
    try
    {
        return GetIpStaticDataTypeMap().find(rIpStaticDataType)->second;
    }
    catch (const std::out_of_range& e)
    {
        throw NuTo::Exception("[NuTo::IpData::IpStaticDataTypeToString] Enum undefined or not implemented.");
    }
}

NuTo::IpData::eIpStaticDataType NuTo::IpData::IpStaticDataTypeToEnum(const std::string& rIpStaticDataType)
{
    std::string uppercase = boost::to_upper_copy(rIpStaticDataType);

    for (auto entry : GetIpStaticDataTypeMap())
        if (entry.second == uppercase)
            return entry.first;

    throw NuTo::Exception("[NuTo::Interpolation::IpStaticDataTypeToEnum] IpStaticDataType " + rIpStaticDataType +
                          " has no enum equivalent or is not implemented.");
}
