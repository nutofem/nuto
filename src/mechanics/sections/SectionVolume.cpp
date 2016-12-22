// $Id$

#include "mechanics/MechanicsException.h"
#include "mechanics/sections/SectionEnum.h"
#include "mechanics/sections/SectionVolume.h"

// constructor
NuTo::SectionVolume::SectionVolume()
{
}

// get section type
NuTo::eSectionType NuTo::SectionVolume::GetType() const
{
    return eSectionType::VOLUME;
}

// info routine
void NuTo::SectionVolume::Info(unsigned short rVerboseLevel) const
{
    this->SectionBase::Info(rVerboseLevel);
}

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::SectionVolume)
#endif
