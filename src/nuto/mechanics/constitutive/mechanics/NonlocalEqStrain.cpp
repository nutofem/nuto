// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/mechanics/NonlocalEqStrain.h"

NuTo::NonlocalEqStrain::NonlocalEqStrain() : ConstitutiveInputBase::ConstitutiveInputBase(), FullVector<double,1>()
{
}


#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::NonlocalEqStrain::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::NonlocalEqStrain::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::NonlocalEqStrain::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::NonlocalEqStrain::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::NonlocalEqStrain::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::NonlocalEqStrain::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::NonlocalEqStrain::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize NonlocalEqStrain" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveInputBase);
    ar & boost::serialization::make_nvp ("NonlocalEqStrainEigen",boost::serialization::base_object< FullVector<double,1> > ( *this ));
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize NonlocalEqStrain" << std::endl;
#endif
}
#endif // ENABLE_SERIALIZATION

