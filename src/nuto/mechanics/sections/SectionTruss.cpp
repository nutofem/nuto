// $Id$

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/sections/SectionTruss.h"
#include "nuto/mechanics/sections/SectionEnum.h"
#include <assert.h>
#include <cmath>

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/vector.hpp>
#endif

// constructor
NuTo::SectionTruss::SectionTruss()
{
    mArea = 0.0;
    mAreaParameters.resize(0);
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
    if (mAreaParameters.size() == 0)
    {
        return 1.;
//        throw NuTo::MechanicsException("[NuTo::SectionTruss::GetAreaFactor] Call NuTo::SectionTruss::SetAreaParameters first!") ;
    }

    double xWeakSpot = mAreaParameters[0];
    double lWeakSpot = mAreaParameters[1];
    double alpha     = mAreaParameters[2];
    //double exponent  = mAreaParameters[3];

    double xStartWeakSpot = xWeakSpot - lWeakSpot / 2;
    double xEndWeakSpot   = xWeakSpot + lWeakSpot / 2;

    double areaFactor = 1.;

    if (xStartWeakSpot < rXCoordinate and rXCoordinate < xEndWeakSpot)
//        areaFactor = 1 - alpha* std::pow(1.-2*std::abs((rXCoordinate-xWeakSpot)/lWeakSpot) ,exponent);
        areaFactor = 1- alpha/2*(1+std::cos(2*M_PI*(rXCoordinate-xWeakSpot) / lWeakSpot));


    assert(areaFactor > 0.);

    return areaFactor;
}


//! @brief ... set the area parameters
//! @param rAreaParameter ... area parameters [0]-xWeakSpot [1]-lWeakSpot [2]-alpha [3]-exponent
void NuTo::SectionTruss::SetAreaParameters(double* rAreaParameters)
{
    for (int i = 0; i < 4; i++) mAreaParameters.push_back(rAreaParameters[i]);
}

// get section type
NuTo::eSectionType NuTo::SectionTruss::GetType() const
{
    return eSectionType::TRUSS;
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

#ifdef ENABLE_SERIALIZATION
template<class Archive>
void NuTo::SectionTruss::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize SectionTruss" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(SectionBase);
    ar & BOOST_SERIALIZATION_NVP(mArea);
    ar & BOOST_SERIALIZATION_NVP(mAreaParameters);
#ifdef DEBUG_SERIALIZATION
    std::cout << "end serialize SectionTruss" << std::endl;
#endif
}


BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::SectionTruss)
#endif // ENABLE_SERIALIZATION
