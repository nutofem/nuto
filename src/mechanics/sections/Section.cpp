#include "base/Exception.h"
#include "mechanics/sections/Section.h"
#include <ostream>

using namespace NuTo;

Section::~Section()
{
}


double Section::GetArea(double) const
{
    throw Exception(__PRETTY_FUNCTION__, "Section type has no cross-section area.");
}


double Section::GetThickness() const
{
    throw Exception(__PRETTY_FUNCTION__, "Section type has no thickness.");
}


double Section::GetCircumference() const
{
    throw Exception(__PRETTY_FUNCTION__, "Section type has no circumference.");
}


bool Section::IsPlaneStrain() const
{
    throw Exception(__PRETTY_FUNCTION__, "Section is not a plane section.");
}


