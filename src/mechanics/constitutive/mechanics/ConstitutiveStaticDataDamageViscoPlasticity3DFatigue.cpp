// $Id: ConstitutiveStaticDataDamageViscoPlasticity3D.cpp 612 2014-01-16 17:04:23Z vkindrac $
// ConstitutiveStaticDataDamageViscoPlasticity3D.cpp

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "mechanics/constitutive/mechanics/ConstitutiveStaticDataDamageViscoPlasticity3DFatigue.h"
#include "mechanics/MechanicsException.h"

//! @brief constructor
NuTo::ConstitutiveStaticDataDamageViscoPlasticity3DFatigue::ConstitutiveStaticDataDamageViscoPlasticity3DFatigue() : ConstitutiveStaticDataDamageViscoPlasticity3D()
{
	mOmegaComprFatigue = 0.;
	mKappaInelasticFatigue = 0.;
	mPrevHardeningFatigue = 0.;

	mEpsilonPFatigue[0] = 0.;
	mEpsilonPFatigue[1] = 0.;
	mEpsilonPFatigue[2] = 0.;
	mEpsilonPFatigue[3] = 0.;
	mEpsilonPFatigue[4] = 0.;
	mEpsilonPFatigue[5] = 0.;

	mEpsilonVpFatigue[0] = 0.;
	mEpsilonVpFatigue[1] = 0.;
	mEpsilonVpFatigue[2] = 0.;
	mEpsilonVpFatigue[3] = 0.;
	mEpsilonVpFatigue[4] = 0.;
	mEpsilonVpFatigue[5] = 0.;

	mPrevStrainFatigue[0] = 0.;
	mPrevStrainFatigue[1] = 0.;
	mPrevStrainFatigue[2] = 0.;
	mPrevStrainFatigue[3] = 0.;
	mPrevStrainFatigue[4] = 0.;
	mPrevStrainFatigue[5] = 0.;

	mPrevStressFatigue[0] = 0.;
	mPrevStressFatigue[1] = 0.;
	mPrevStressFatigue[2] = 0.;
	mPrevStressFatigue[3] = 0.;
	mPrevStressFatigue[4] = 0.;
	mPrevStressFatigue[5] = 0.;
}

NuTo::ConstitutiveStaticDataDamageViscoPlasticity3DFatigue& NuTo::ConstitutiveStaticDataDamageViscoPlasticity3DFatigue::operator= (NuTo::ConstitutiveStaticDataDamageViscoPlasticity3DFatigue const& rOther)
{
    NuTo::ConstitutiveStaticDataDamageViscoPlasticity3D::operator= (rOther);
    mOmegaComprFatigue = rOther.mOmegaComprFatigue;
    mKappaInelasticFatigue = rOther.mKappaInelasticFatigue;
    mPrevHardeningFatigue = rOther.mPrevHardeningFatigue;

    mEpsilonPFatigue[0] = rOther.mEpsilonPFatigue[0];
    mEpsilonPFatigue[1] = rOther.mEpsilonPFatigue[1];
    mEpsilonPFatigue[2] = rOther.mEpsilonPFatigue[2];
    mEpsilonPFatigue[3] = rOther.mEpsilonPFatigue[3];
    mEpsilonPFatigue[4] = rOther.mEpsilonPFatigue[4];
    mEpsilonPFatigue[5] = rOther.mEpsilonPFatigue[5];

    mEpsilonVpFatigue[0] = rOther.mEpsilonVpFatigue[0];
    mEpsilonVpFatigue[1] = rOther.mEpsilonVpFatigue[1];
    mEpsilonVpFatigue[2] = rOther.mEpsilonVpFatigue[2];
    mEpsilonVpFatigue[3] = rOther.mEpsilonVpFatigue[3];
    mEpsilonVpFatigue[4] = rOther.mEpsilonVpFatigue[4];
    mEpsilonVpFatigue[5] = rOther.mEpsilonVpFatigue[5];

    mPrevStrainFatigue[0] = rOther.mPrevStrainFatigue[0];
    mPrevStrainFatigue[1] = rOther.mPrevStrainFatigue[1];
    mPrevStrainFatigue[2] = rOther.mPrevStrainFatigue[2];
    mPrevStrainFatigue[3] = rOther.mPrevStrainFatigue[3];
    mPrevStrainFatigue[4] = rOther.mPrevStrainFatigue[4];
    mPrevStrainFatigue[5] = rOther.mPrevStrainFatigue[5];

    mPrevStressFatigue[0] = rOther.mPrevStressFatigue[0];
    mPrevStressFatigue[1] = rOther.mPrevStressFatigue[1];
    mPrevStressFatigue[2] = rOther.mPrevStressFatigue[2];
    mPrevStressFatigue[3] = rOther.mPrevStressFatigue[3];
    mPrevStressFatigue[4] = rOther.mPrevStressFatigue[4];
    mPrevStressFatigue[5] = rOther.mPrevStressFatigue[5];
    return *this;
}


//! @brief check, if the static data is compatible with a given element and a given constitutive model
bool NuTo::ConstitutiveStaticDataDamageViscoPlasticity3DFatigue::CheckConstitutiveCompatibility(NuTo::Constitutive::eConstitutiveType rConstitutiveType, NuTo::Element::eElementType rElementType)const
{
	if (rConstitutiveType==NuTo::Constitutive::DAMAGE_VISCO_PLASTICITY_HARDENING_ENGINEERING_STRESS ||
			rConstitutiveType==NuTo::Constitutive::DAMAGE_VISCO_PLASTICITY_ENGINEERING_STRESS)
	{
		if (rElementType==NuTo::Element::BRICK8N || rElementType==NuTo::Element::TETRAHEDRON4N || rElementType==NuTo::Element::TETRAHEDRON10N)
			return true;
		else
			return false;
	}
	else
		return false;
}

//!@ brief reinterpret as damage viscoplasticity 3D static data
NuTo::ConstitutiveStaticDataDamageViscoPlasticity3DFatigue* NuTo::ConstitutiveStaticDataDamageViscoPlasticity3DFatigue::AsDamageViscoPlasticity3DFatigue()
{
    return this;
}

//!@ brief reinterpret as damage viscoplasticity 3D static data
const NuTo::ConstitutiveStaticDataDamageViscoPlasticity3DFatigue* NuTo::ConstitutiveStaticDataDamageViscoPlasticity3DFatigue::AsDamageViscoPlasticity3DFatigue()const
{
    return this;
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::ConstitutiveStaticDataDamageViscoPlasticity3DFatigue::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataDamageViscoPlasticity3DFatigue::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataDamageViscoPlasticity3DFatigue::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataDamageViscoPlasticity3DFatigue::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataDamageViscoPlasticity3DFatigue::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataDamageViscoPlasticity3DFatigue::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstitutiveStaticDataDamageViscoPlasticity3DFatigue::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstitutiveStaticDataDamageViscoPlasticity3D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveStaticDataDamageViscoPlasticity3D)
       & BOOST_SERIALIZATION_NVP(mOmegaComprFatigue)
       & BOOST_SERIALIZATION_NVP(mKappaInelasticFatigue)
       & BOOST_SERIALIZATION_NVP(mEpsilonPFatigue)
       & BOOST_SERIALIZATION_NVP(mEpsilonVpFatigue)
       & BOOST_SERIALIZATION_NVP(mPrevStrainFatigue)
       & BOOST_SERIALIZATION_NVP(mPrevStressFatigue)
       & BOOST_SERIALIZATION_NVP(mPrevHardeningFatigue);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstitutiveStaticDataDamageViscoPlasticity3DFatigue" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstitutiveStaticDataDamageViscoPlasticity3DFatigue)
#endif // ENABLE_SERIALIZATION
