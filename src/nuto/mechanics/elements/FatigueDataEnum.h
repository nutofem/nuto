// $Id: FatigueDataEnum.h 2014-12-16 vkindrac $
#ifndef FATIGUEDATAENUM_H_
#define FATIGUEDATAENUM_H_

namespace NuTo
{
namespace FatigueData
{

//! @brief covers all fatigue data

enum eFatigueDataType
{
    SAVE_STATIC_DATA,          //!< save static data in case if the cycle jump should be repeated
    RESTORE_STATIC_DATA,       //!< cycle jump should be repeated, getting the previous state
    EXTRAPOLATE_STATIC_DATA    //!< extrapolate static data to make a jump
};

}
}
#endif /* FATIGUEDATAENUM_H_ */
