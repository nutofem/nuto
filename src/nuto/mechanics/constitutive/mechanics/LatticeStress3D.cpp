// $Id: $

#include <iostream>

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/mechanics/LatticeStress3D.h"

NuTo::LatticeStress3D::LatticeStress3D()
{
    this->mLatticeStress[0] = 0.0;
    this->mLatticeStress[1] = 0.0;
    this->mLatticeStress[2] = 0.0;
}

unsigned int NuTo::LatticeStress3D::GetNumberOfComponents() const
{
    return 3;
}

const double* NuTo::LatticeStress3D::GetData() const
{
    return this->mLatticeStress;
}

void NuTo::LatticeStress3D::Info(unsigned short rVerboseLevel) const
{
    std::cout << "    components of the lattice traction vector: "
              << this->mLatticeStress[0] << ", "
              << this->mLatticeStress[1] << ", "
              << this->mLatticeStress[2] << std::endl;
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::LatticeStress3D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::LatticeStress3D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::LatticeStress3D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::LatticeStress3D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::LatticeStress3D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::LatticeStress3D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::LatticeStress3D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize LatticeStress3D" << std::endl;
#endif
   ar & BOOST_SERIALIZATION_NVP(mLatticeStress);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize LatticeStress3D" << std::endl;
#endif
}
#endif // ENABLE_SERIALIZATION
