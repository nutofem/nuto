// $Id$

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/sections/SectionTruss.h"
#include <assert.h>
#include <cmath>

// constructor
NuTo::SectionTruss::SectionTruss()
{
    this->mArea = 0.0;
    this->mAreaParameters = nullptr;
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

//! @brief ... calculates a x-dependent area factor based on the mAreaParameters data
//! @param rXCoordinate ... x coordinate
//! @return ... x-dependent area factor
double NuTo::SectionTruss::GetAreaFactor(double rXCoordinate) const
{
    if (mAreaParameters == nullptr)
    {
        throw NuTo::MechanicsException("[NuTo::SectionTruss::GetAreaFactor] Call NuTo::SectionTruss::SetAreaParameters first!") ;
    }

    double xWeakSpot = mAreaParameters[0];
    double lWeakSpot = mAreaParameters[1];
    double alpha     = mAreaParameters[2];
    double exponent  = mAreaParameters[3];

    double xStartWeakSpot = xWeakSpot - lWeakSpot / 2;
    double xEndWeakSpot   = xWeakSpot + lWeakSpot / 2;

    double areaFactor = 1.;

    if (xStartWeakSpot < rXCoordinate and rXCoordinate < xEndWeakSpot)
        return 1 - alpha* std::pow(1.-std::abs((rXCoordinate-xWeakSpot)/lWeakSpot*2.) ,exponent);

    assert(areaFactor > 0.);

    return areaFactor;
}


//! @brief ... set the area parameters
//! @param rAreaParameter ... area parameters [0]-xWeakSpot [1]-lWeakSpot [2]-alpha [3]-exponent
void NuTo::SectionTruss::SetAreaParameters(double* rAreaParameters)
{
    mAreaParameters = rAreaParameters;
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

//! @brief ... cast the base pointer to an SectionTruss, otherwise throws an exception
NuTo::SectionTruss* NuTo::SectionTruss::AsSectionTruss()
{
    return this;
}
//! @brief ... cast the base pointer to an SectionTruss, otherwise throws an exception
const NuTo::SectionTruss* NuTo::SectionTruss::AsSectionTruss() const
{
    return this;
}
