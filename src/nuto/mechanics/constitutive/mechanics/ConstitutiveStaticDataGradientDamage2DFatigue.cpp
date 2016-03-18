
#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataGradientDamage2DFatigue.h"
#include "nuto/mechanics/MechanicsException.h"

//! @brief constructor
NuTo::ConstitutiveStaticDataGradientDamage2DFatigue::ConstitutiveStaticDataGradientDamage2DFatigue() : ConstitutiveStaticDataPrevEngineeringStressStrain2DPlaneStrain()
{
    mKappaFatigue = 0;
    mOmegaFatigue = 0;

    mKappaExtrapolated = 0;
    mOmegaExtrapolated = 0.;

    mKappaDeltaImplicit = 0.;
    mOmegaDeltaImplicit = 0.;

    mKappa = 0.;
    mOmega = 0;

	mPrevStrainFatigue[0] = 0.;
	mPrevStrainFatigue[1] = 0.;
	mPrevStrainFatigue[2] = 0.;

	mPrevStressFatigue[0] = 0.;
	mPrevStressFatigue[1] = 0.;
	mPrevStressFatigue[2] = 0.;

    mPrevNonlocalEqStrain.SetValue(0,0.);
    mPrevNonlocalEqStrainFatigue.SetValue(0,0.);

}

NuTo::ConstitutiveStaticDataGradientDamage2DFatigue& NuTo::ConstitutiveStaticDataGradientDamage2DFatigue::operator= (NuTo::ConstitutiveStaticDataGradientDamage2DFatigue const& rOther)
{
    NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain2DPlaneStrain::operator= (rOther);
    mKappaFatigue = rOther.mKappaFatigue;
    mOmegaFatigue = rOther.mOmegaFatigue;
    mKappaExtrapolated = rOther.mKappaExtrapolated;
    mOmegaExtrapolated = rOther.mOmegaExtrapolated;
    mKappaDeltaImplicit = rOther.mKappaDeltaImplicit;
    mOmegaDeltaImplicit = rOther.mOmegaDeltaImplicit;
    mKappa        = rOther.mKappa;
    mOmega        = rOther.mOmega;

    mPrevStrainFatigue[0] = rOther.mPrevStrainFatigue[0];
    mPrevStrainFatigue[1] = rOther.mPrevStrainFatigue[1];
    mPrevStrainFatigue[2] = rOther.mPrevStrainFatigue[2];

    mPrevStressFatigue[0] = rOther.mPrevStressFatigue[0];
    mPrevStressFatigue[1] = rOther.mPrevStressFatigue[1];
    mPrevStressFatigue[2] = rOther.mPrevStressFatigue[2];

    mPrevNonlocalEqStrain = rOther.mPrevNonlocalEqStrain;
    mPrevNonlocalEqStrainFatigue = rOther.mPrevNonlocalEqStrainFatigue;

    return *this;
}

NuTo::ConstitutiveStaticDataGradientDamage2DFatigue* NuTo::ConstitutiveStaticDataGradientDamage2DFatigue::Clone()const
{
	return new ConstitutiveStaticDataGradientDamage2DFatigue(*this);
}


#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::ConstitutiveStaticDataGradientDamage2DFatigue::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataGradientDamage2DFatigue::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataGradientDamage2DFatigue::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataGradientDamage2DFatigue::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataGradientDamage2DFatigue::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataGradientDamage2DFatigue::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstitutiveStaticDataGradientDamage2DFatigue::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstitutiveStaticDataGradientDamage2DFatigue" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveStaticDataPrevEngineeringStressStrain2DPlaneStrain)
       & BOOST_SERIALIZATION_NVP(mKappaFatigue)
	   & BOOST_SERIALIZATION_NVP(mOmegaFatigue)
       & BOOST_SERIALIZATION_NVP(mKappaExtrapolated)
	   & BOOST_SERIALIZATION_NVP(mOmegaExtrapolated)
       & BOOST_SERIALIZATION_NVP(mKappaDeltaImplicit)
	   & BOOST_SERIALIZATION_NVP(mOmegaDeltaImplicit)
       & BOOST_SERIALIZATION_NVP(mKappa)
	   & BOOST_SERIALIZATION_NVP(mOmega)
       & BOOST_SERIALIZATION_NVP(mPrevStrainFatigue)
       & BOOST_SERIALIZATION_NVP(mPrevStressFatigue)
	   & BOOST_SERIALIZATION_NVP(mPrevNonlocalEqStrain)
	   & BOOST_SERIALIZATION_NVP(mPrevNonlocalEqStrainFatigue);
;
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstitutiveStaticDataGradientDamage2DFatigue" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstitutiveStaticDataGradientDamage2DFatigue)
#endif // ENABLE_SERIALIZATION

//!@ brief reinterpret as nonlocal damage2D static data
NuTo::ConstitutiveStaticDataGradientDamage2DFatigue* NuTo::ConstitutiveStaticDataGradientDamage2DFatigue::AsGradientDamage2DFatigue()
{
    return this;
}

//!@ brief reinterpret as nonlocal damage2d static data
const NuTo::ConstitutiveStaticDataGradientDamage2DFatigue* NuTo::ConstitutiveStaticDataGradientDamage2DFatigue::AsGradientDamage2DFatigue()const
{
    return this;
}

//! @brief check, if the static data is compatible with a given element and a given constitutive model
bool NuTo::ConstitutiveStaticDataGradientDamage2DFatigue::CheckConstitutiveCompatibility(NuTo::Constitutive::eConstitutiveType rConstitutiveType, NuTo::Element::eElementType rElementType)const
{
	if (rConstitutiveType==NuTo::Constitutive::GRADIENT_DAMAGE_ENGINEERING_STRESS ||
			rConstitutiveType==NuTo::Constitutive::GRADIENT_DAMAGE_ENGINEERING_STRESS_FATIGUE)
	{
//		if (rElementType==NuTo::Element::ELEMENT2D)
		if (rElementType==NuTo::Element::PLANE2D6N) {
			// delete cout
			std::cout << " CSD:: Place 22 " << std::endl;
			return true;}
		else
			{
			// delete cout
			std::cout << " CSD:: Place 222 " << std::endl;
			return false;
			}
	}
	else
		return false;
}

