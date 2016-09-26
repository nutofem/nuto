// $Id$ 
#pragma once

#include <map>
#include <string>

namespace NuTo
{
namespace IpData
{
enum class eIpDataType
{
    NOIPDATA,                //!< no additional ip data
    STATICDATA,              //!< static data
    STATICDATANONLOCAL,      //!< nonlocal and static data
    MULTISCALE               //!< multiscale - a full structure on the fine scale whose average values are used
};

const std::map<eIpDataType, std::string> GetIpDataTypeMap();

//! @brief covers all ip data (not only static data) that is dependent on the current iteration state
//! @brief this is mainly used in Get routines for visualization purposes
enum class eIpStaticDataType
{
    BOND_STRESS,                    //!< bond stress
    DAMAGE,                         //!< isotropic damage variable
    ELASTIC_ENERGY,                 //!< elastic energy
    ENGINEERING_PLASTIC_STRAIN,     //!< plastic strain
    ENGINEERING_STRAIN,             //!< engineering strain
    ENGINEERING_STRESS,             //!< engineering stress
    EXTRAPOLATION_ERROR,            //!< for implicit / explicit time integration schemes
    HEAT_FLUX,                      //!< heat flux
    INTERNAL_ENERGY,                //!< internal (elastic + inelastic) energy
    LATTICE_STRAIN,                 //!< lattice strain
    LATTICE_STRESS,                 //!< lattice stress
    LATTICE_PLASTIC_STRAIN,         //!< lattice plastic strain
    LOCAL_EQ_STRAIN,                //!< local equivalent strain
    SHRINKAGE_STRAIN,               //!< shrinkage strains
    THERMAL_STRAIN,                 //!< thermal strains
    SLIP,                           //!< slip, i.e. relative displacement
    TOTAL_INELASTIC_EQ_STRAIN       //!< total inelastic equivalent strain
};



const std::map<eIpStaticDataType, std::string> GetIpStaticDataTypeMap();
std::string IpDataTypeToString(const eIpDataType& rIpDataType);
std::string IpStaticDataTypeToString(const eIpStaticDataType& rIpStaticDataType);
eIpDataType IpDataTypeToEnum(const std::string& rIpDataType);
eIpStaticDataType IpStaticDataTypeToEnum(const std::string& rIpStaticDataType);





}
}

