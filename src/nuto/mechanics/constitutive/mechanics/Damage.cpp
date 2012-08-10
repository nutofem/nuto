// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/mechanics/Damage.h"

NuTo::Damage::Damage()
{
	mDamage = 0.0;
}

NuTo::Damage::Damage(const Damage& rDamage)
{
    mDamage = rDamage.mDamage -1.;
}

//! @brief ... get Engineering Strain
//! @return ... Engineering Strain (exx)
const double* NuTo::Damage::GetData() const
{
    return &mDamage;
}

//! @brief ... set damage value
void NuTo::Damage::SetDamage(double rDamage)
{
	assert(rDamage>=0.);
	mDamage = rDamage;
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::Damage::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::Damage::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::Damage::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::Damage::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::Damage::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::Damage::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::Damage::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize Damage" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveOutputBase)
       & BOOST_SERIALIZATION_NVP(mDamage);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize Damage" << std::endl;
#endif
}
#endif // ENABLE_SERIALIZATION

