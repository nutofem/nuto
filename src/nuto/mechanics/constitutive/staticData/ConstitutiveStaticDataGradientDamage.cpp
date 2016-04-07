#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/staticData/ConstitutiveStaticDataGradientDamage.h"

//! @brief constructor
NuTo::ConstitutiveStaticDataGradientDamage::ConstitutiveStaticDataGradientDamage() :
        NuTo::ConstitutiveStaticDataBase::ConstitutiveStaticDataBase()
{
    mKappa = 0;
}

NuTo::ConstitutiveStaticDataGradientDamage* NuTo::ConstitutiveStaticDataGradientDamage::Clone() const
{
    return new ConstitutiveStaticDataGradientDamage(*this);
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::ConstitutiveStaticDataGradientDamage::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataGradientDamage::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataGradientDamage::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataGradientDamage::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataGradientDamage::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataGradientDamage::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstitutiveStaticDataGradientDamage::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstitutiveStaticDataGradientDamage" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveStaticDataBase)
       & BOOST_SERIALIZATION_NVP(mKappa);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstitutiveStaticDataGradientDamage" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstitutiveStaticDataGradientDamage)
#endif // ENABLE_SERIALIZATION

//!@ brief reinterpret as nonlocal damage1d static data
NuTo::ConstitutiveStaticDataGradientDamage* NuTo::ConstitutiveStaticDataGradientDamage::AsGradientDamage()
{
    return this;
}

//!@ brief reinterpret as nonlocal damage2d static data
const NuTo::ConstitutiveStaticDataGradientDamage* NuTo::ConstitutiveStaticDataGradientDamage::AsGradientDamage() const
{
    return this;
}

//! @brief check, if the static data is compatible with a given element and a given constitutive model
bool NuTo::ConstitutiveStaticDataGradientDamage::CheckConstitutiveCompatibility(NuTo::Constitutive::eConstitutiveType rConstitutiveType, NuTo::Element::eElementType rElementType) const
{
    if (rConstitutiveType == NuTo::Constitutive::GRADIENT_DAMAGE_ENGINEERING_STRESS)
    {
        if (rElementType == NuTo::Element::CONTINUUMELEMENT)
            return true;
        if (rElementType == NuTo::Element::CONTINUUMBOUNDARYELEMENT)
            return true;
    }
    return false;
}

