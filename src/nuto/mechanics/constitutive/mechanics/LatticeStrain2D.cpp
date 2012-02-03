// $Id: LatticeStrain2D.cpp 316 2010-09-28 19:40:50Z unger3 $

#include <iostream>

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/mechanics/LatticeStrain2D.h"

NuTo::LatticeStrain2D::LatticeStrain2D()
{
    this->mLatticeStrain[0] = 0.0;
    this->mLatticeStrain[1] = 0.0;
}

unsigned int NuTo::LatticeStrain2D::GetNumberOfComponents() const
{
    return 2;
}

const double* NuTo::LatticeStrain2D::GetData() const
{
    return this->mLatticeStrain;
}

void NuTo::LatticeStrain2D::Info(unsigned short rVerboseLevel) const
{
    std::cout << "    components of the lattice strain vector: "
              << this->mLatticeStrain[0] << ", "
              << this->mLatticeStrain[1] << std::endl;
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::LatticeStrain2D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::LatticeStrain2D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::LatticeStrain2D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::LatticeStrain2D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::LatticeStrain2D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::LatticeStrain2D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::LatticeStrain2D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize LatticeStrain2D" << std::endl;
#endif
   ar & BOOST_SERIALIZATION_NVP(mLatticeStrain);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize LatticeStrain2D" << std::endl;
#endif
}
#endif // ENABLE_SERIALIZATION
