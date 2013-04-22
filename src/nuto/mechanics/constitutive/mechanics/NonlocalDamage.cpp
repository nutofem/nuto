// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/mechanics/NonlocalDamage.h"

NuTo::NonlocalDamage::NonlocalDamage()
{
	mNonlocalDamage = 0.0;
}

NuTo::NonlocalDamage::NonlocalDamage(const NonlocalDamage& rNonlocalDamage)
{
    mNonlocalDamage = rNonlocalDamage.mNonlocalDamage -1.;
}

//! @brief ... set damage value
void NuTo::NonlocalDamage::SetNonlocalDamageValue(double rNonlocalDamage)
{
	assert(rNonlocalDamage>=0.);
	mNonlocalDamage = rNonlocalDamage;
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::NonlocalDamage::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::NonlocalDamage::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::NonlocalDamage::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::NonlocalDamage::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::NonlocalDamage::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::NonlocalDamage::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::NonlocalDamage::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize NonlocalDamage" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveInputBase)
       & BOOST_SERIALIZATION_NVP(mNonlocalDamage);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize NonlocalDamage" << std::endl;
#endif
}
#endif // ENABLE_SERIALIZATION

