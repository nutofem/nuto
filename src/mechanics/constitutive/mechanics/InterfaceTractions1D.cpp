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

#include "mechanics/constitutive/mechanics/InterfaceTractions1D.h"

NuTo::InterfaceTractions1D::InterfaceTractions1D()
{
    this->mInterfaceTractions = 0.0;
}

unsigned int NuTo::InterfaceTractions1D::GetNumberOfComponents() const
{
    return 1;
}

const double* NuTo::InterfaceTractions1D::GetData() const
{
    return &mInterfaceTractions;
}


void NuTo::InterfaceTractions1D::Info(unsigned short rVerboseLevel) const
{
    std::cout << "    components of the interface traction vector: " << this->mInterfaceTractions << std::endl;
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::InterfaceTractions1D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::InterfaceTractions1D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::InterfaceTractions1D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::InterfaceTractions1D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::InterfaceTractions1D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::InterfaceTractions1D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::InterfaceTractions1D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize InterfaceTractions1D" << std::endl;
#endif
   ar & BOOST_SERIALIZATION_NVP(mInterfaceTractions);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize InterfaceTractions1D" << std::endl;
#endif
}
#endif // ENABLE_SERIALIZATION
