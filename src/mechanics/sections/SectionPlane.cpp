#include <iostream>
#include "mechanics/MechanicsException.h"
#include "mechanics/sections/SectionPlane.h"

using namespace NuTo;

SectionPlane::SectionPlane(double thickness, bool isPlaneStrain)
    : mThickness(thickness)
    , mIsPlaneStrain(isPlaneStrain)
{
}


std::shared_ptr<SectionPlane> SectionPlane::Create(double thickness, bool isPlaneStrain)
{
    return std::shared_ptr<SectionPlane>(new SectionPlane(thickness, isPlaneStrain));
}


double SectionPlane::GetThickness() const
{
    return mThickness;
}


void SectionPlane::Info(std::ostream& out) const
{
    out << "    Plane section with thickness: " << mThickness << "\n";
    if (mIsPlaneStrain)
        out << "    Section type is plane strain.\n";
    else
        out << "    Section type is plane stress.\n";
}


bool SectionPlane::IsPlaneStrain() const
{
    return mIsPlaneStrain;
}

