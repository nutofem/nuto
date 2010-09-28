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

#include "nuto/mechanics/constitutive/mechanics/InterfaceTractions2D.h"

NuTo::InterfaceTractions2D::InterfaceTractions2D()
{
    this->mInterfaceTractions[0] = 0.0;
    this->mInterfaceTractions[1] = 0.0;
}

unsigned int NuTo::InterfaceTractions2D::GetNumberOfComponents() const
{
    return 2;
}

const double* NuTo::InterfaceTractions2D::GetData() const
{
    return this->mInterfaceTractions;
}

void NuTo::InterfaceTractions2D::Info(unsigned short rVerboseLevel) const
{
    std::cout << "    components of the interface traction vector: "
              << this->mInterfaceTractions[0] << ", "
              << this->mInterfaceTractions[1] << std::endl;
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::InterfaceTractions2D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::InterfaceTractions2D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::InterfaceTractions2D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::InterfaceTractions2D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::InterfaceTractions2D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::InterfaceTractions2D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::InterfaceTractions2D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize InterfaceTractions2D" << std::endl;
#endif
   ar & BOOST_SERIALIZATION_NVP(mInterfaceTractions);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize InterfaceTractions2D" << std::endl;
#endif
}
#endif // ENABLE_SERIALIZATION
