#include <iostream>
#include "mechanics/sections/SectionFibreMatrixBond.h"

using namespace NuTo;


SectionFibreMatrixBond::SectionFibreMatrixBond(double circumference):
mCircumference(circumference)
{
}


std::shared_ptr<SectionFibreMatrixBond> SectionFibreMatrixBond::Create(double circumference)
{
    return std::shared_ptr<SectionFibreMatrixBond>(new SectionFibreMatrixBond(circumference));
}


double SectionFibreMatrixBond::GetCircumference() const
{
    return mCircumference;
}


void SectionFibreMatrixBond::Info(std::ostream& out) const
{
    out << "    Fibre matrix bond section with circumference: " << mCircumference << "\n";
}

