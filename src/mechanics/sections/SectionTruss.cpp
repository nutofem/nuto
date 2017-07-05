#include <iostream>
#include "mechanics/MechanicsException.h"
#include "mechanics/sections/SectionTruss.h"

using namespace NuTo;

SectionTruss::SectionTruss(double area)
    : mArea(area)
{
}

std::shared_ptr<SectionTruss> SectionTruss::Create(double area)
{
    return std::shared_ptr<SectionTruss>(new SectionTruss(area));
}


double SectionTruss::GetArea(double) const
{
    return mArea;
}


void SectionTruss::Info(std::ostream& out) const
{
    out << "    Truss section with area " << mArea << "\n";
}
