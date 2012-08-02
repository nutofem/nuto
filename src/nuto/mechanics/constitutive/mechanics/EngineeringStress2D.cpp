// $Id$
#include <iostream>

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/mechanics/EngineeringStress2D.h"

// constructor
NuTo::EngineeringStress2D::EngineeringStress2D(): ConstitutiveOutputBase::ConstitutiveOutputBase()
{
    this->mEngineeringStress[0] = 0.0;
    this->mEngineeringStress[1] = 0.0;
    this->mEngineeringStress[2] = 0.0;
}

// number of components
unsigned int NuTo::EngineeringStress2D::GetNumberOfComponents() const
{
    return 3;
}

// get Engineering stress
const double* NuTo::EngineeringStress2D::GetData() const
{
    return this->mEngineeringStress;
}

// info routine
void NuTo::EngineeringStress2D::Info(unsigned short rVerboseLevel) const
{
    std::cout << "    components of Engineering stress tensor (vector notation): " << this->mEngineeringStress[0] << ", "
              << this->mEngineeringStress[1] << ", " << this->mEngineeringStress[2]  << std::endl;
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::EngineeringStress2D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::EngineeringStress2D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::EngineeringStress2D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::EngineeringStress2D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::EngineeringStress2D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::EngineeringStress2D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::EngineeringStress2D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize EngineeringStress2D" << std::endl;
#endif
    ar &  BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveOutputBase)
       & BOOST_SERIALIZATION_NVP(mEngineeringStress);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize EngineeringStress2D" << std::endl;
#endif
}
#endif // ENABLE_SERIALIZATION
