
#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "mechanics/constitutive/mechanics/ConstitutiveStaticDataGradientDamage1DFatigue.h"
#include "mechanics/MechanicsException.h"

//! @brief constructor
NuTo::ConstitutiveStaticDataGradientDamage1DFatigue::ConstitutiveStaticDataGradientDamage1DFatigue() : ConstitutiveStaticDataPrevEngineeringStressStrain1D()
{
    mKappaFatigue = 0;
    mOmegaFatigue = 0;
    mKappa = 0.;
    mOmega = 0;
    mPrevNonlocalEqStrain.SetValue(0,0.);
}

NuTo::ConstitutiveStaticDataGradientDamage1DFatigue& NuTo::ConstitutiveStaticDataGradientDamage1DFatigue::operator= (NuTo::ConstitutiveStaticDataGradientDamage1DFatigue const& rOther)
{
    NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain1D::operator= (rOther);
    mKappaFatigue = rOther.mKappaFatigue;
    mOmegaFatigue = rOther.mOmegaFatigue;
    mKappa        = rOther.mKappa;
    mOmega        = rOther.mOmega;
    mPrevNonlocalEqStrain = rOther.mPrevNonlocalEqStrain;
    return *this;
}

NuTo::ConstitutiveStaticDataGradientDamage1DFatigue* NuTo::ConstitutiveStaticDataGradientDamage1DFatigue::Clone()const
{
	return new ConstitutiveStaticDataGradientDamage1DFatigue(*this);
}


#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::ConstitutiveStaticDataGradientDamage1DFatigue::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataGradientDamage1DFatigue::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataGradientDamage1DFatigue::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataGradientDamage1DFatigue::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataGradientDamage1DFatigue::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataGradientDamage1DFatigue::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstitutiveStaticDataGradientDamage1DFatigue::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstitutiveStaticDataGradientDamage1DFatigue" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveStaticDataPrevEngineeringStressStrain1D)
       & BOOST_SERIALIZATION_NVP(mKappaFatigue)
	   & BOOST_SERIALIZATION_NVP(mOmegaFatigue)
       & BOOST_SERIALIZATION_NVP(mKappa)
	   & BOOST_SERIALIZATION_NVP(mOmega)
	   & BOOST_SERIALIZATION_NVP(mPrevNonlocalEqStrain);
;
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstitutiveStaticDataGradientDamage1DFatigue" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstitutiveStaticDataGradientDamage1DFatigue)
#endif // ENABLE_SERIALIZATION

//!@ brief reinterpret as nonlocal damage1d static data
NuTo::ConstitutiveStaticDataGradientDamage1DFatigue* NuTo::ConstitutiveStaticDataGradientDamage1DFatigue::AsGradientDamage1DFatigue()
{
    return this;
}

//!@ brief reinterpret as nonlocal damage2d static data
const NuTo::ConstitutiveStaticDataGradientDamage1DFatigue* NuTo::ConstitutiveStaticDataGradientDamage1DFatigue::AsGradientDamage1DFatigue()const
{
    return this;
}

//! @brief check, if the static data is compatible with a given element and a given constitutive model
bool NuTo::ConstitutiveStaticDataGradientDamage1DFatigue::CheckConstitutiveCompatibility(NuTo::Constitutive::eConstitutiveType rConstitutiveType, NuTo::Element::eElementType rElementType)const
{
	if (rConstitutiveType==NuTo::Constitutive::GRADIENT_DAMAGE_ENGINEERING_STRESS ||
			rConstitutiveType==NuTo::Constitutive::GRADIENT_DAMAGE_ENGINEERING_STRESS_FATIGUE)
	{
		if (rElementType==NuTo::Element::TRUSS1D3N)
			return true;
		else
			return false;
	}
	else
		return false;
}

