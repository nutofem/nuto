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

#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataDamageViscoPlasticity3D.h"
#include "nuto/mechanics/MechanicsException.h"

//! @brief constructor
NuTo::ConstitutiveStaticDataDamageViscoPlasticity3D::ConstitutiveStaticDataDamageViscoPlasticity3D() : ConstitutiveStaticDataPrevEngineeringStressStrain3D()
{
	mOmegaP = 0.;
	mKappaP = 0.;
	mVP     = 0.;

	mEpsilonP[0] = 0.;
	mEpsilonP[1] = 0.;
	mEpsilonP[2] = 0.;
	mEpsilonP[3] = 0.;
	mEpsilonP[4] = 0.;
	mEpsilonP[5] = 0.;
}

NuTo::ConstitutiveStaticDataDamageViscoPlasticity3D& NuTo::ConstitutiveStaticDataDamageViscoPlasticity3D::operator= (NuTo::ConstitutiveStaticDataDamageViscoPlasticity3D const& rOther)
{
    NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain3D::operator= (rOther);
    mOmegaP = rOther.mOmegaP;
    mKappaP= rOther.mKappaP;
    mVP = rOther.mVP;
    mEpsilonP[0] = rOther.mEpsilonP[0];
    mEpsilonP[1] = rOther.mEpsilonP[1];
    mEpsilonP[2] = rOther.mEpsilonP[2];
    mEpsilonP[3] = rOther.mEpsilonP[3];
    mEpsilonP[4] = rOther.mEpsilonP[4];
    mEpsilonP[5] = rOther.mEpsilonP[5];
    return *this;
}


//! @brief check, if the static data is compatible with a given element and a given constitutive model
bool NuTo::ConstitutiveStaticDataDamageViscoPlasticity3D::CheckConstitutiveCompatibility(NuTo::Constitutive::eConstitutiveType rConstitutiveType, NuTo::Element::eElementType rElementType)const
{
	if (rConstitutiveType==NuTo::Constitutive::NONLOCAL_DAMAGE_PLASTICITY_ENGINEERING_STRESS)
	{
		if (rElementType==NuTo::Element::BRICK8N || rElementType==NuTo::Element::TETRAHEDRON4N || rElementType==NuTo::Element::TETRAHEDRON10N)
			return true;
		else
			return false;
	}
	else
		return false;
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::ConstitutiveStaticDataDamageViscoPlasticity3D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataDamageViscoPlasticity3D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataDamageViscoPlasticity3D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataDamageViscoPlasticity3D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataDamageViscoPlasticity3D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataDamageViscoPlasticity3D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstitutiveStaticDataDamageViscoPlasticity3D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstitutiveStaticDataDamageViscoPlasticity3D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveStaticDataPrevEngineeringStressStrain3D)
       & BOOST_SERIALIZATION_NVP(mOmega)
       & BOOST_SERIALIZATION_NVP(mKappa)
       & BOOST_SERIALIZATION_NVP(mEpsilonP);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstitutiveStaticDataDamageViscoPlasticity3D" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstitutiveStaticDataDamageViscoPlasticity3D)
#endif // ENABLE_SERIALIZATION
