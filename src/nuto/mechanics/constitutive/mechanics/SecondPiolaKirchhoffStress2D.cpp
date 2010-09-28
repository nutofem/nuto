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

#include "nuto/mechanics/constitutive/mechanics/SecondPiolaKirchhoffStress2D.h"

// constructor
NuTo::SecondPiolaKirchhoffStress2D::SecondPiolaKirchhoffStress2D()
{
    this->mSecondPiolaKirchhoffStress[0] = 0.0;
    this->mSecondPiolaKirchhoffStress[1] = 0.0;
    this->mSecondPiolaKirchhoffStress[2] = 0.0;
}

// number of components
unsigned int NuTo::SecondPiolaKirchhoffStress2D::GetNumberOfComponents() const
{
    return 3;
}

// get second Piola-Kirchhoff stress
const double* NuTo::SecondPiolaKirchhoffStress2D::GetData() const
{
    return this->mSecondPiolaKirchhoffStress;
}

// info routine
void NuTo::SecondPiolaKirchhoffStress2D::Info(unsigned short rVerboseLevel) const
{
    std::cout << "    components of second Piola-Kirchhoff stress tensor (vector notation): "
              << this->mSecondPiolaKirchhoffStress[0] << ", "
              << this->mSecondPiolaKirchhoffStress[1] << ", "
              << this->mSecondPiolaKirchhoffStress[2] << std::endl;
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::SecondPiolaKirchhoffStress2D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::SecondPiolaKirchhoffStress2D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::SecondPiolaKirchhoffStress2D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::SecondPiolaKirchhoffStress2D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::SecondPiolaKirchhoffStress2D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::SecondPiolaKirchhoffStress2D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::SecondPiolaKirchhoffStress2D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize SecondPiolaKirchhoffStress2D" << std::endl;
#endif
   ar & BOOST_SERIALIZATION_NVP(mSecondPiolaKirchhoffStress);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize SecondPiolaKirchhoffStress2D" << std::endl;
#endif
}
#endif // ENABLE_SERIALIZATION
