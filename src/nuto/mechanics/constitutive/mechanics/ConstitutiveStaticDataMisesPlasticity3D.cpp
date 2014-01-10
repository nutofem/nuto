// $Id: ConstitutiveStaticDataBase.cpp 102 2009-11-11 10:47:23Z eckardt4 $
#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataMisesPlasticity3D.h"
#include "nuto/mechanics/MechanicsException.h"

//! @brief constructor
NuTo::ConstitutiveStaticDataMisesPlasticity3D::ConstitutiveStaticDataMisesPlasticity3D() : ConstitutiveStaticDataBase()
{
	mEpsilonPEq = 0.;

	mEpsilonP[0] = 0.;
	mEpsilonP[1] = 0.;
	mEpsilonP[2] = 0.;
	mEpsilonP[3] = 0.;
	mEpsilonP[4] = 0.;
	mEpsilonP[5] = 0.;

    mSigmaB[0] = 0.;
    mSigmaB[1] = 0.;
    mSigmaB[2] = 0.;
    mSigmaB[3] = 0.;
    mSigmaB[4] = 0.;
    mSigmaB[5] = 0.;
}

//! @brief assignment operator
NuTo::ConstitutiveStaticDataMisesPlasticity3D& NuTo::ConstitutiveStaticDataMisesPlasticity3D::operator= (ConstitutiveStaticDataMisesPlasticity3D const& rOther)
{
    NuTo::ConstitutiveStaticDataBase::operator= (rOther);
	mEpsilonPEq = rOther.mEpsilonPEq;

	mEpsilonP[0] = rOther.mEpsilonP[0] ;
	mEpsilonP[1] = rOther.mEpsilonP[1] ;
	mEpsilonP[2] = rOther.mEpsilonP[2];
	mEpsilonP[3] = rOther.mEpsilonP[3];
	mEpsilonP[4] = rOther.mEpsilonP[4];
	mEpsilonP[5] = rOther.mEpsilonP[5];

    mSigmaB[0] = rOther.mSigmaB[0];
    mSigmaB[1] = rOther.mSigmaB[1] ;
    mSigmaB[2] = rOther.mSigmaB[2];
    mSigmaB[3] = rOther.mSigmaB[3];
    mSigmaB[4] = rOther.mSigmaB[4];
    mSigmaB[5] = rOther.mSigmaB[5];

    return (*this);
}

//! @brief check, if the static data is compatible with a given element and a given constitutive model
bool NuTo::ConstitutiveStaticDataMisesPlasticity3D::CheckConstitutiveCompatibility(NuTo::Constitutive::eConstitutiveType rConstitutiveType, NuTo::Element::eElementType rElementType)const
{
	if (rConstitutiveType==NuTo::Constitutive::MISES_PLASTICITY_ENGINEERING_STRESS)
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
template void NuTo::ConstitutiveStaticDataMisesPlasticity3D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataMisesPlasticity3D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataMisesPlasticity3D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataMisesPlasticity3D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataMisesPlasticity3D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataMisesPlasticity3D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstitutiveStaticDataMisesPlasticity3D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstitutiveStaticDataMisesPlasticity3D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveStaticDataBase)
       & BOOST_SERIALIZATION_NVP(mEpsilonPEq)
       & BOOST_SERIALIZATION_NVP(mEpsilonP)
       & BOOST_SERIALIZATION_NVP(mSigmaB);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstitutiveStaticDataMisesPlasticity3D" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstitutiveStaticDataMisesPlasticity3D)
#endif // ENABLE_SERIALIZATION
