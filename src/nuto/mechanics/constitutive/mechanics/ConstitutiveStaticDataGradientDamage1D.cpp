
#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataGradientDamage1D.h"
#include "nuto/mechanics/MechanicsException.h"

//! @brief constructor
NuTo::ConstitutiveStaticDataGradientDamage1D::ConstitutiveStaticDataGradientDamage1D()
   : NuTo::ConstitutiveStaticDataBase::ConstitutiveStaticDataBase()
{
    mKappa = 0;
}

NuTo::ConstitutiveStaticDataGradientDamage1D& NuTo::ConstitutiveStaticDataGradientDamage1D::operator= (NuTo::ConstitutiveStaticDataGradientDamage1D const& rOther)
{
    NuTo::ConstitutiveStaticDataBase::operator= (rOther);
    mKappa = rOther.mKappa;
    return *this;
}

NuTo::ConstitutiveStaticDataGradientDamage1D* NuTo::ConstitutiveStaticDataGradientDamage1D::Clone()const
{
	return new ConstitutiveStaticDataGradientDamage1D(*this);
}


#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::ConstitutiveStaticDataGradientDamage1D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataGradientDamage1D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataGradientDamage1D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataGradientDamage1D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataGradientDamage1D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataGradientDamage1D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstitutiveStaticDataGradientDamage1D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstitutiveStaticDataGradientDamage1D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveStaticDataBase)
       & BOOST_SERIALIZATION_NVP(mKappa);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstitutiveStaticDataGradientDamage1D" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstitutiveStaticDataGradientDamage1D)
#endif // ENABLE_SERIALIZATION

//!@ brief reinterpret as nonlocal damage1d static data
NuTo::ConstitutiveStaticDataGradientDamage1D* NuTo::ConstitutiveStaticDataGradientDamage1D::AsGradientDamage1D()
{
    return this;
}

//!@ brief reinterpret as nonlocal damage2d static data
const NuTo::ConstitutiveStaticDataGradientDamage1D* NuTo::ConstitutiveStaticDataGradientDamage1D::AsGradientDamage1D()const
{
    return this;
}

//! @brief check, if the static data is compatible with a given element and a given constitutive model
bool NuTo::ConstitutiveStaticDataGradientDamage1D::CheckConstitutiveCompatibility(NuTo::Constitutive::eConstitutiveType rConstitutiveType, NuTo::Element::eElementType rElementType)const
{
	if (rConstitutiveType==NuTo::Constitutive::GRADIENT_DAMAGE_ENGINEERING_STRESS)
	{
		if (rElementType==NuTo::Element::TRUSS1D3N)
			return true;
		else
			return false;
	}
	else
		return false;
}

