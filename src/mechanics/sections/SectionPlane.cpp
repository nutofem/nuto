#include <iostream>
#include "mechanics/sections/SectionPlane.h"

using namespace NuTo;

SectionPlane::SectionPlane(double thickness, bool isPlaneStrain)
    : mThickness(thickness)
    , mIsPlaneStrain(isPlaneStrain)
	, mIsAxiSymmetric(false)
{
}

SectionPlane::SectionPlane()
    : mThickness(1.)
    , mIsPlaneStrain(true)
	, mIsAxiSymmetric(true)
{
}

std::shared_ptr<SectionPlane> SectionPlane::Create(double thickness, bool isPlaneStrain)
{
	return std::shared_ptr<SectionPlane>(new SectionPlane(thickness, isPlaneStrain));
}

std::shared_ptr<SectionPlane> SectionPlane::CreateAxiSymmetric()
{
	return std::shared_ptr<SectionPlane>(new SectionPlane());
}

double SectionPlane::GetThickness() const
{
    return mThickness;
}


void SectionPlane::Info(std::ostream& out) const
{
    out << "    Plane section with thickness: " << mThickness << "\n";
    if (mIsAxiSymmetric) {
    	out << "    Section type is axisymmetric.\n";
	} else {
	    if (mIsPlaneStrain)
	        out << "    Section type is plane strain.\n";
	    else
	        out << "    Section type is plane stress.\n";
	}
}

bool SectionPlane::IsPlaneStrain() const
{
    return mIsPlaneStrain;
}

bool SectionPlane::IsAxiSymmetric() const
{
    return mIsAxiSymmetric;
}
