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

#include "nuto/mechanics/constitutive/mechanics/EngineeringStress1D.h"

// constructor
NuTo::EngineeringStress1D::EngineeringStress1D(): ConstitutiveOutputBase::ConstitutiveOutputBase()
{
	(*this)[0] = 0.0;
}

// number of components
unsigned int NuTo::EngineeringStress1D::GetNumberOfComponents() const
{
    return 1;
}

// get Engineering stress
const double* NuTo::EngineeringStress1D::GetData() const
{
    return data();
}

// info routine
void NuTo::EngineeringStress1D::Info(unsigned short rVerboseLevel) const
{
    std::cout << "    components of Engineering stress tensor (vector notation): " << (*this)[0] << std::endl;
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::EngineeringStress1D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::EngineeringStress1D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::EngineeringStress1D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::EngineeringStress1D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::EngineeringStress1D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::EngineeringStress1D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::EngineeringStress1D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize EngineeringStress1D" << std::endl;
#endif
    ar &  BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveOutputBase);
    ar & boost::serialization::make_nvp ("EngineeringStress1DEigen",boost::serialization::base_object< FullVectorFixed<double,1> > ( *this ));
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize EngineeringStress1D" << std::endl;
#endif
}
#endif // ENABLE_SERIALIZATION
