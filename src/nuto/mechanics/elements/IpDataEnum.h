// $Id$ 
#ifndef IPDATAENUM_H_
#define IPDATAENUM_H_

#include <map>
#include <boost/algorithm/string.hpp>
#include "nuto/mechanics/MechanicsException.h"

namespace NuTo
{
namespace IpData
{
enum eIpDataType
{
    NOIPDATA,                //!< no additional ip data
    STATICDATA,              //!< static data
    STATICDATANONLOCAL,      //!< nonlocal and static data
    STATICDATAWEIGHTCOORDINATES2D, //!< static data, local weight and coordinates stored at the integration point level
    STATICDATAWEIGHTCOORDINATES3D, //!< static data, local weight and coordinates stored at the integration point level
    MULTISCALE               //!< multiscale - a full structure on the fine scale whose average values are used
};
static inline std::map<eIpDataType, std::string> GetIpDataTypeMap()
{
    std::map<eIpDataType, std::string> shapeTypeMap;
    shapeTypeMap[eIpDataType::NOIPDATA]                      = "NOIPDATA";
    shapeTypeMap[eIpDataType::STATICDATA]                    = "STATICDATA";
    shapeTypeMap[eIpDataType::STATICDATANONLOCAL]            = "STATICDATANONLOCAL";
    shapeTypeMap[eIpDataType::STATICDATAWEIGHTCOORDINATES2D] = "STATICDATAWEIGHTCOORDINATES2D";
    shapeTypeMap[eIpDataType::STATICDATAWEIGHTCOORDINATES3D] = "STATICDATAWEIGHTCOORDINATES3D";
    shapeTypeMap[eIpDataType::MULTISCALE]                    = "MULTISCALE";
    return shapeTypeMap;
}

//! @brief covers all ip data (not only static data) that is dependent on the current iteration state
//! @brief this is mainly used in Get routines for visualization purposes
enum eIpStaticDataType
{
    LATTICE_STRAIN,            //!< lattice strain
    LATTICE_STRESS,            //!< lattice stress
    LATTICE_PLASTIC_STRAIN,    //!< lattice plastic strain
    ENGINEERING_STRAIN,        //!< engineering strain
    ENGINEERING_STRESS,        //!< engineering stress
    LOCAL_EQ_STRAIN,           //!< local equivalent strain
    DAMAGE,                    //!< isotropic damage variable
    ENGINEERING_PLASTIC_STRAIN,//!> plastic strain
    ELASTIC_ENERGY,            //!> elastic energy
    INTERNAL_ENERGY,           //!> internal (elastic + inelastic) energy
    HEAT_FLUX                  //!> heat flux
};
static inline std::map<eIpStaticDataType, std::string> GetIpStaticDataTypeMap()
{
    std::map<eIpStaticDataType, std::string> shapeTypeMap;
    shapeTypeMap[eIpStaticDataType::LATTICE_STRAIN]             = "LATTICE_STRAIN";
    shapeTypeMap[eIpStaticDataType::LATTICE_STRESS]             = "LATTICE_STRESS";
    shapeTypeMap[eIpStaticDataType::LATTICE_PLASTIC_STRAIN]     = "LATTICE_PLASTIC_STRAIN";
    shapeTypeMap[eIpStaticDataType::ENGINEERING_STRAIN]         = "ENGINEERING_STRAIN";
    shapeTypeMap[eIpStaticDataType::ENGINEERING_STRESS]         = "ENGINEERING_STRESS";
    shapeTypeMap[eIpStaticDataType::DAMAGE]                     = "DAMAGE";
    shapeTypeMap[eIpStaticDataType::ENGINEERING_PLASTIC_STRAIN] = "ENGINEERING_PLASTIC_STRAIN";
    shapeTypeMap[eIpStaticDataType::ELASTIC_ENERGY]             = "ELASTIC_ENERGY";
    shapeTypeMap[eIpStaticDataType::INTERNAL_ENERGY]            = "INTERNAL_ENERGY";
    shapeTypeMap[eIpStaticDataType::HEAT_FLUX]                  = "HEAT_FLUX";
    return shapeTypeMap;
}





static inline std::string IpDataTypeToString(const eIpDataType& rIpDataType)
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

static inline std::string IpStaticDataTypeToString(const eIpStaticDataType& rIpStaticDataType)
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

static inline eIpDataType IpDataTypeToEnum(const std::string& rIpDataType)
{
    std::string uppercase = boost::to_upper_copy(rIpDataType);

    for(auto entry : GetIpDataTypeMap())
        if (entry.second == uppercase)
            return entry.first;

    throw NuTo::MechanicsException("[NuTo::Interpolation::IpDataTypeToEnum] IpDataType " + rIpDataType + " has no enum equivalent or is not implemented.");
}

static inline eIpStaticDataType IpStaticDataTypeToEnum(const std::string& rIpStaticDataType)
{
    std::string uppercase = boost::to_upper_copy(rIpStaticDataType);

    for(auto entry : GetIpStaticDataTypeMap())
        if (entry.second == uppercase)
            return entry.first;

    throw NuTo::MechanicsException("[NuTo::Interpolation::IpStaticDataTypeToEnum] IpStaticDataType " + rIpStaticDataType + " has no enum equivalent or is not implemented.");
}





}
}
#endif /* IPDATAENUM_H_ */
