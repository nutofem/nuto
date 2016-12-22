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

#include "mechanics/constitutive/mechanics/InterfaceTractions3D.h"

NuTo::InterfaceTractions3D::InterfaceTractions3D()
{
    this->mInterfaceTractions[0] = 0.0;
    this->mInterfaceTractions[1] = 0.0;
    this->mInterfaceTractions[2] = 0.0;
}

unsigned int NuTo::InterfaceTractions3D::GetNumberOfComponents() const
{
    return 3;
}

const double* NuTo::InterfaceTractions3D::GetData() const
{
    return this->mInterfaceTractions;
}

void NuTo::InterfaceTractions3D::Info(unsigned short rVerboseLevel) const
{
    std::cout << "    components of the interface traction vector: "
              << this->mInterfaceTractions[0] << ", "
              << this->mInterfaceTractions[1] << ", "
              << this->mInterfaceTractions[2] << std::endl;
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::InterfaceTractions3D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::InterfaceTractions3D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::InterfaceTractions3D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::InterfaceTractions3D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::InterfaceTractions3D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::InterfaceTractions3D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::InterfaceTractions3D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize InterfaceTractions3D" << std::endl;
#endif
   ar & BOOST_SERIALIZATION_NVP(mInterfaceTractions);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize InterfaceTractions3D" << std::endl;
#endif
}
#endif // ENABLE_SERIALIZATION
