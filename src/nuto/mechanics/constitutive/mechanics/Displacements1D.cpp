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

#include "nuto/mechanics/constitutive/mechanics/Displacements1D.h"

// constructor
NuTo::Displacements1D::Displacements1D()
{
    this->mDisplacements = 0.0;
}

// get displacements
const double* NuTo::Displacements1D::GetData() const
{
    return &mDisplacements;
}

// get displacements
double* NuTo::Displacements1D::GetData()
{
    return &mDisplacements;
}

// info routine
void NuTo::Displacements1D::Info(unsigned short rVerboseLevel) const
{
    std::cout << "    displacement: " << this->mDisplacements << std::endl;
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::Displacements1D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::Displacements1D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::Displacements1D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::Displacements1D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::Displacements1D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::Displacements1D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::Displacements1D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize Displacements1D" << std::endl;
#endif
   ar & BOOST_SERIALIZATION_NVP(mDisplacements);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize Displacements1D" << std::endl;
#endif
}
#endif // ENABLE_SERIALIZATION
