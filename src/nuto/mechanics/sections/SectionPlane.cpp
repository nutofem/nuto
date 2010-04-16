// $Id$

#include "nuto/mechanics/MechanicsException.h"
#include <nuto/mechanics/sections/SectionPlane.h>

// constructor
NuTo::SectionPlane::SectionPlane(eSectionType rSectionType)
{
    this->mThickness = 0.0;
    if (rSectionType==PLANE_STRAIN || rSectionType==PLANE_STRESS)
    	this->mSectionType = rSectionType;
    else
    {
    	throw NuTo::MechanicsException("[NuTo::SectionPlane::SectionPlane] section for plane elements is either plane strain or plane stress.") ;
    }
}

// set section thickness
void NuTo::SectionPlane::SetThickness(double rThickness)
{
    if (rThickness <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::SectionPlane::SetThickness] section thickness must be greater than 0.") ;
    }
    this->mThickness = rThickness;
}

// get section thickness
double NuTo::SectionPlane::GetThickness() const
{
    return this->mThickness;
}

// get section type
NuTo::SectionBase::eSectionType NuTo::SectionPlane::GetType() const
{
    return mSectionType;
}

// info routine
void NuTo::SectionPlane::Info(unsigned short rVerboseLevel) const
{
    this->SectionBase::Info(rVerboseLevel);
    std::cout << "    section type: 2D" << std::endl;
    std::cout << "    section thickness: " << this->mThickness << std::endl;
    if (mSectionType==PLANE_STRAIN)
        std::cout << "    section stress state: plane strain" << std::endl;
    else
    {
        if (mSectionType==PLANE_STRESS)
            std::cout << "    section stress state: plane stress" << std::endl;
    }
}
