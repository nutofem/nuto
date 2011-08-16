// $Id$ 
// ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain.cpp
// created May 6, 2010 by Joerg F. Unger

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain.h"
#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataPrevEngineeringStressStrain2DPlaneStrain.h"
#include "nuto/mechanics/MechanicsException.h"

//! @brief constructor
NuTo::ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain::ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain()
   : NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain2DPlaneStrain::ConstitutiveStaticDataPrevEngineeringStressStrain2DPlaneStrain()
{
    mKappa = 0.;

	mEpsilonP[0] = 0.;
	mEpsilonP[1] = 0.;
	mEpsilonP[2] = 0.;
	mEpsilonP[3] = 0.;

	mTmpdEpsilonPdEpsilon[0] = 0.;
	mTmpdEpsilonPdEpsilon[1] = 0.;
	mTmpdEpsilonPdEpsilon[2] = 0.;
	mTmpdEpsilonPdEpsilon[3] = 0.;
	mTmpdEpsilonPdEpsilon[4] = 0.;
	mTmpdEpsilonPdEpsilon[5] = 0.;
	mTmpdEpsilonPdEpsilon[6] = 0.;
	mTmpdEpsilonPdEpsilon[7] = 0.;
	mTmpdEpsilonPdEpsilon[8] = 0.;
	mTmpdEpsilonPdEpsilon[9] = 0.;
	mTmpdEpsilonPdEpsilon[10] = 0.;
	mTmpdEpsilonPdEpsilon[11] = 0.;
	mTmpdEpsilonPdEpsilon[12] = 0.;
	mTmpdEpsilonPdEpsilon[13] = 0.;
	mTmpdEpsilonPdEpsilon[14] = 0.;
	mTmpdEpsilonPdEpsilon[15] = 0.;

	mTmpKappa = 0.;

	mTmpEpsilonP[0] = 0.;
	mTmpEpsilonP[1] = 0.;
	mTmpEpsilonP[2] = 0.;
	mTmpEpsilonP[3] = 0.;

	mTmpLeq = 0.;

	mTmpdLeqdEpsilon[0] = 0.;
	mTmpdLeqdEpsilon[1] = 0.;
	mTmpdLeqdEpsilon[2] = 0.;
	mTmpdLeqdEpsilon[3] = 0.;
}

NuTo::ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain& NuTo::ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain::operator= (NuTo::ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain const& rOther)
{
    NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain2DPlaneStrain::operator= (rOther);
    mKappa = rOther.mKappa;
    mEpsilonP[0] = rOther.mEpsilonP[0];
    mEpsilonP[1] = rOther.mEpsilonP[1];
    mEpsilonP[2] = rOther.mEpsilonP[2];
    mEpsilonP[3] = rOther.mEpsilonP[3];
	mTmpdEpsilonPdEpsilon[0] = rOther.mTmpdEpsilonPdEpsilon[0];
	mTmpdEpsilonPdEpsilon[1] = rOther.mTmpdEpsilonPdEpsilon[1];
	mTmpdEpsilonPdEpsilon[2] = rOther.mTmpdEpsilonPdEpsilon[2];
	mTmpdEpsilonPdEpsilon[3] = rOther.mTmpdEpsilonPdEpsilon[3];
	mTmpdEpsilonPdEpsilon[4] = rOther.mTmpdEpsilonPdEpsilon[4];
	mTmpdEpsilonPdEpsilon[5] = rOther.mTmpdEpsilonPdEpsilon[5];
	mTmpdEpsilonPdEpsilon[6] = rOther.mTmpdEpsilonPdEpsilon[6];
	mTmpdEpsilonPdEpsilon[7] = rOther.mTmpdEpsilonPdEpsilon[7];
	mTmpdEpsilonPdEpsilon[8] = rOther.mTmpdEpsilonPdEpsilon[8];
	mTmpdEpsilonPdEpsilon[9] = rOther.mTmpdEpsilonPdEpsilon[9];
	mTmpdEpsilonPdEpsilon[10] = rOther.mTmpdEpsilonPdEpsilon[10];
	mTmpdEpsilonPdEpsilon[11] = rOther.mTmpdEpsilonPdEpsilon[11];
	mTmpdEpsilonPdEpsilon[12] = rOther.mTmpdEpsilonPdEpsilon[12];
	mTmpdEpsilonPdEpsilon[13] = rOther.mTmpdEpsilonPdEpsilon[13];
	mTmpdEpsilonPdEpsilon[14] = rOther.mTmpdEpsilonPdEpsilon[14];
	mTmpdEpsilonPdEpsilon[15] = rOther.mTmpdEpsilonPdEpsilon[15];
    mTmpKappa = rOther.mTmpKappa;
    mTmpEpsilonP[0] = rOther.mTmpEpsilonP[0];
    mTmpEpsilonP[1] = rOther.mTmpEpsilonP[1];
    mTmpEpsilonP[2] = rOther.mTmpEpsilonP[2];
    mTmpEpsilonP[3] = rOther.mTmpEpsilonP[3];
    mTmpLeq = rOther.mTmpLeq;
    mTmpdLeqdEpsilon[0] = rOther.mTmpdLeqdEpsilon[0];
    mTmpdLeqdEpsilon[1] = rOther.mTmpdLeqdEpsilon[1];
    mTmpdLeqdEpsilon[2] = rOther.mTmpdLeqdEpsilon[2];
    mTmpdLeqdEpsilon[3] = rOther.mTmpdLeqdEpsilon[3];
    return *this;
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveStaticDataPrevEngineeringStressStrain2DPlaneStrain)
       & BOOST_SERIALIZATION_NVP(mKappa)
       & BOOST_SERIALIZATION_NVP(mEpsilonP)
       & BOOST_SERIALIZATION_NVP(mTmpdEpsilonPdEpsilon)
       & BOOST_SERIALIZATION_NVP(mTmpKappa)
       & BOOST_SERIALIZATION_NVP(mTmpEpsilonP)
       & BOOST_SERIALIZATION_NVP(mTmpLeq)
       & BOOST_SERIALIZATION_NVP(mTmpdLeqdEpsilon);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain)
#endif // ENABLE_SERIALIZATION

//!@ brief reinterpret as nonlocal damage2d static data
NuTo::ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain* NuTo::ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain::AsNonlocalDamagePlasticity2DPlaneStrain()
{
    return this;
}

//!@ brief reinterpret as nonlocal damage2d static data
const NuTo::ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain* NuTo::ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain::AsNonlocalDamagePlasticity2DPlaneStrain()const
{
    return this;
}

//! @brief check, if the static data is compatible with a given element and a given constitutive model
bool NuTo::ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain::CheckConstitutiveCompatibility(NuTo::Constitutive::eConstitutiveType rConstitutiveType, NuTo::Element::eElementType rElementType)const
{
	if (rConstitutiveType==NuTo::Constitutive::NONLOCAL_DAMAGE_PLASTICITY)
	{
		if (rElementType==NuTo::Element::PLANE2D3N || rElementType==NuTo::Element::PLANE2D4N || rElementType==NuTo::Element::PLANE2D6N)
			return true;
		else
			return false;
	}
	else
		return false;
}

