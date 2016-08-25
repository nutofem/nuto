#include "IpDataEnum.h"

#include <boost/algorithm/string.hpp>
#include "nuto/mechanics/MechanicsException.h"

const std::map<NuTo::IpData::eIpDataType, std::string> NuTo::IpData::GetIpDataTypeMap()
{
    const std::map<eIpDataType, std::string> shapeTypeMap =
       {{eIpDataType::NOIPDATA,             "NOIPDATA"},
        {eIpDataType::STATICDATA,           "STATICDATA"},
        {eIpDataType::STATICDATANONLOCAL,   "STATICDATANONLOCAL"},
        {eIpDataType::MULTISCALE,           "MULTISCALE"}};

    return shapeTypeMap;
}

const std::map<NuTo::IpData::eIpStaticDataType, std::string> NuTo::IpData::GetIpStaticDataTypeMap()
{
    const std::map<eIpStaticDataType, std::string> shapeTypeMap =
       {{eIpStaticDataType::DAMAGE,     "DAMAGE"},
        {eIpStaticDataType::ELASTIC_ENERGY, "ELASTIC_ENERGY"},
        {eIpStaticDataType::ENGINEERING_PLASTIC_STRAIN, "ENGINEERING_PLASTIC_STRAIN"},
        {eIpStaticDataType::ENGINEERING_STRAIN, "ENGINEERING_STRAIN"},
        {eIpStaticDataType::ENGINEERING_STRESS, "ENGINEERING_STRESS"},
        {eIpStaticDataType::EXTRAPOLATION_ERROR, "EXTRAPOLATION_ERROR"},
        {eIpStaticDataType::HEAT_FLUX, "HEAT_FLUX"},
        {eIpStaticDataType::INTERNAL_ENERGY, "INTERNAL_ENERGY"},
        {eIpStaticDataType::LATTICE_PLASTIC_STRAIN, "LATTICE_PLASTIC_STRAIN"},
        {eIpStaticDataType::LATTICE_STRAIN, "LATTICE_STRAIN"},
        {eIpStaticDataType::LATTICE_STRESS, "LATTICE_STRESS"},
        {eIpStaticDataType::SHRINKAGE_STRAIN, "SHRINKAGE_STRAIN"},
        {eIpStaticDataType::THERMAL_STRAIN, "THERMAL_STRAIN"},
        {eIpStaticDataType::TOTAL_INELASTIC_EQ_STRAIN, "TOTAL_INELASTIC_EQUIVALENT_STRAIN"}};
    return shapeTypeMap;
}

std::string NuTo::IpData::IpDataTypeToString(const NuTo::IpData::eIpDataType &rIpDataType)
{
    try
    {
        return GetIpDataTypeMap().find(rIpDataType)->second;
    }
    catch (const std::out_of_range& e)
    {
        throw NuTo::MechanicsException("[NuTo::IpData::IpDataTypeToString] Enum undefined or not implemented.");
    }
}

std::string NuTo::IpData::IpStaticDataTypeToString(const NuTo::IpData::eIpStaticDataType &rIpStaticDataType)
{
    try
    {
        return GetIpStaticDataTypeMap().find(rIpStaticDataType)->second;
    }
    catch (const std::out_of_range& e)
    {
        throw NuTo::MechanicsException("[NuTo::IpData::IpStaticDataTypeToString] Enum undefined or not implemented.");
    }
}

NuTo::IpData::eIpDataType NuTo::IpData::IpDataTypeToEnum(const std::string &rIpDataType)
{
    std::string uppercase = boost::to_upper_copy(rIpDataType);

    for(auto entry : GetIpDataTypeMap())
        if (entry.second == uppercase)
            return entry.first;

    throw NuTo::MechanicsException("[NuTo::Interpolation::IpDataTypeToEnum] IpDataType " + rIpDataType + " has no enum equivalent or is not implemented.");
}

NuTo::IpData::eIpStaticDataType NuTo::IpData::IpStaticDataTypeToEnum(const std::string &rIpStaticDataType)
{
    std::string uppercase = boost::to_upper_copy(rIpStaticDataType);

    for(auto entry : GetIpStaticDataTypeMap())
        if (entry.second == uppercase)
            return entry.first;

    throw NuTo::MechanicsException("[NuTo::Interpolation::IpStaticDataTypeToEnum] IpStaticDataType " + rIpStaticDataType + " has no enum equivalent or is not implemented.");
}
