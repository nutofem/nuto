#include <iostream>
#include <cassert>
#include <cmath>
#include "mechanics/sections/SectionVariableTruss.h"

using namespace NuTo;

SectionVariableTruss::SectionVariableTruss(double area, double location, double extent, double reductionRatio)
    : SectionTruss(area)
    , mWeakSpotLocation(location)
    , mWeakSpotExtent(extent)
    , mAlpha(reductionRatio)
{
}

std::shared_ptr<SectionVariableTruss> SectionVariableTruss::Create(double area, double location, double extent,
                                                                   double reductionRatio)
{
    return std::shared_ptr<SectionVariableTruss>(new SectionVariableTruss(area, location, extent, reductionRatio));
}


double SectionVariableTruss::GetArea(double coordinate) const
{
    double start = mWeakSpotLocation - mWeakSpotExtent / 2;
    double end = mWeakSpotLocation + mWeakSpotExtent / 2;

    if (start < coordinate and coordinate < end)
    {
        double areaFactor =
                1 - mAlpha / 2 * (1 + std::cos(2 * M_PI * (coordinate - mWeakSpotLocation) / mWeakSpotExtent));
        assert(areaFactor > 0.);
        return areaFactor * SectionTruss::GetArea(coordinate);
    }
    else
        return SectionTruss::GetArea(coordinate);
}


void SectionVariableTruss::Info() const
{
    std::cout << "    Variable truss section." << std::endl;
}
