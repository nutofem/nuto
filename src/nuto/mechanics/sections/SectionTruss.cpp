// $Id$

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/sections/SectionTruss.h"

// constructor
NuTo::SectionTruss::SectionTruss()
{
    this->mArea = 0.0;
}

// set cross-section area
void NuTo::SectionTruss::SetArea(double rArea)
{
    if (rArea <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::SectionTruss::SetArea] cross-section area must greater than 0.") ;
    }
    this->mArea = rArea;
}

// get cross-section area
double NuTo::SectionTruss::GetArea() const
{
    return this->mArea;
}

// get section type
NuTo::Section::eSectionType NuTo::SectionTruss::GetType() const
{
    return Section::TRUSS;
}

// info routine
void NuTo::SectionTruss::Info(unsigned short rVerboseLevel) const
{
    this->SectionBase::Info(rVerboseLevel);
    std::cout << "    section type: 1D" << std::endl;
    std::cout << "    cross-section area: " << this->mArea << std::endl;
}
