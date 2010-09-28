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

#include "nuto/mechanics/constitutive/mechanics/SecondPiolaKirchhoffStress1D.h"

// constructor
NuTo::SecondPiolaKirchhoffStress1D::SecondPiolaKirchhoffStress1D()
{
    this->mSecondPiolaKirchhoffStress = 0.0;
}

// number of components
unsigned int NuTo::SecondPiolaKirchhoffStress1D::GetNumberOfComponents() const
{
    return 1;
}

// get second Piola-Kirchhoff stress
const double* NuTo::SecondPiolaKirchhoffStress1D::GetData() const
{
    return &mSecondPiolaKirchhoffStress;
}

// info routine
void NuTo::SecondPiolaKirchhoffStress1D::Info(unsigned short rVerboseLevel) const
{
    std::cout << "    components of second Piola-Kirchhoff stress tensor (vector notation): " << this->mSecondPiolaKirchhoffStress << std::endl;
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::SecondPiolaKirchhoffStress1D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::SecondPiolaKirchhoffStress1D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::SecondPiolaKirchhoffStress1D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::SecondPiolaKirchhoffStress1D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::SecondPiolaKirchhoffStress1D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::SecondPiolaKirchhoffStress1D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::SecondPiolaKirchhoffStress1D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize SecondPiolaKirchhoffStress1D" << std::endl;
#endif
   ar & BOOST_SERIALIZATION_NVP(mSecondPiolaKirchhoffStress);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize SecondPiolaKirchhoffStress1D" << std::endl;
#endif
}
#endif // ENABLE_SERIALIZATION
