// $Id$ 
#ifndef IPDATAENUM_H_
#define IPDATAENUM_H_

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

//! @brief covers all ip data (not only static data) that is dependent on the current iteration state
//! @brief this is mainly used in Get routines for visualization purposes
enum eIpStaticDataType
{
    LATTICE_STRAIN,            //!< lattice strain
    LATTICE_STRESS,            //!< lattice stress
    LATTICE_PLASTIC_STRAIN,    //!< lattice plastic strain
    ENGINEERING_STRAIN,        //!< engineering strain
    ENGINEERING_STRESS,        //!< engineering stress
    DAMAGE,                    //!< isotropic damage variable
    ENGINEERING_PLASTIC_STRAIN,//!> plastic strain
    ELASTIC_ENERGY,            //!> elastic energy
    INTERNAL_ENERGY,           //!> internal (elastic + inelastic) energy
    HEAT_FLUX                  //!> heat flux
};

}
}
#endif /* IPDATAENUM_H_ */
