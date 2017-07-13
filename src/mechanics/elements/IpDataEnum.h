#pragma once

#include <map>
#include <string>

namespace NuTo
{
namespace IpData
{

//! @brief covers all ip data (not only static data) that is dependent on the current iteration state
//! @brief this is mainly used in Get routines for visualization purposes
enum class eIpStaticDataType
{
    BOND_STRESS, //!< bond stress
    DAMAGE, //!< isotropic damage variable
    ELECTRIC_FIELD, //!< electric field calculated from potential
    ELECTRIC_DISPLACEMENT, //!< (di)electric displacment
    ELASTIC_ENERGY, //!< elastic energy
    ENGINEERING_PLASTIC_STRAIN, //!< plastic strain
    ENGINEERING_STRAIN, //!< engineering strain
    ENGINEERING_STRESS, //!< engineering stress
    EXTRAPOLATION_ERROR, //!< for implicit / explicit time integration schemes
    HEAT_FLUX, //!< heat flux
    INTERNAL_ENERGY, //!< internal (elastic + inelastic) energy
    LATTICE_STRAIN, //!< lattice strain
    LATTICE_STRESS, //!< lattice stress
    LATTICE_PLASTIC_STRAIN, //!< lattice plastic strain
    LOCAL_EQ_STRAIN, //!< local equivalent strain
    SHRINKAGE_STRAIN, //!< shrinkage strains
    THERMAL_STRAIN, //!< thermal strains
    SLIP, //!< slip, i.e. relative displacement
    TOTAL_INELASTIC_EQ_STRAIN //!< total inelastic equivalent strain
};


const std::map<eIpStaticDataType, std::string> GetIpStaticDataTypeMap();
std::string IpStaticDataTypeToString(const eIpStaticDataType& rIpStaticDataType);
eIpStaticDataType IpStaticDataTypeToEnum(const std::string& rIpStaticDataType);
}
}
