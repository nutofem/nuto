// $Id$

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/sections/SectionVolume.h"

// constructor
NuTo::SectionVolume::SectionVolume()
{
}

// get section type
NuTo::Section::eSectionType NuTo::SectionVolume::GetType() const
{
    return Section::VOLUME;
}

// info routine
void NuTo::SectionVolume::Info(unsigned short rVerboseLevel) const
{
    this->SectionBase::Info(rVerboseLevel);
}
