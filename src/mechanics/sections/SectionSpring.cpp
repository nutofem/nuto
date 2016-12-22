// $Id$

#include "mechanics/sections/SectionSpring.h"
#include "mechanics/sections/SectionEnum.h"
#include <iostream>

// constructor
NuTo::SectionSpring::SectionSpring()
{
}


// get section type
NuTo::eSectionType NuTo::SectionSpring::GetType() const
{
    return eSectionType::SPRING;
}

// info routine
void NuTo::SectionSpring::Info(unsigned short rVerboseLevel) const
{
    this->SectionBase::Info(rVerboseLevel);
    std::cout << "    section type: 1D" << std::endl;
}

