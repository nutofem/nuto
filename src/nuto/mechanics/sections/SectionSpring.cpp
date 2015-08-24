// $Id$

#include "nuto/mechanics/sections/SectionSpring.h"


// constructor
NuTo::SectionSpring::SectionSpring()
{
}


// get section type
NuTo::Section::eSectionType NuTo::SectionSpring::GetType() const
{
    return Section::SPRING;
}

// info routine
void NuTo::SectionSpring::Info(unsigned short rVerboseLevel) const
{
    this->SectionBase::Info(rVerboseLevel);
    std::cout << "    section type: 1D" << std::endl;
}

