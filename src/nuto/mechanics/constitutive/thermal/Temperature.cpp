// $Id: EngineeringStress1D.cpp 316 2010-09-28 19:40:50Z unger3 $
#include <iostream>

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/thermal/Temperature.h"

// constructor
NuTo::Temperature::Temperature(): ConstitutiveInputBase()
{
    this->mTemperature = 0.0;
}

// get temperature
const double& NuTo::Temperature::GetData() const
{
    return mTemperature;
}

// get temperature
double& NuTo::Temperature::GetData()
{
    return mTemperature;
}

// info routine
void NuTo::Temperature::Info(unsigned short rVerboseLevel) const
{
    std::cout << "    temperature: " << this->mTemperature << std::endl;
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::Temperature::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::Temperature::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::Temperature::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::Temperature::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::Temperature::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::Temperature::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::Temperature::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize Temperature" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveInputBase)
       & BOOST_SERIALIZATION_NVP(mTemperature);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize Temperature" << std::endl;
#endif
}
#endif // ENABLE_SERIALIZATION
