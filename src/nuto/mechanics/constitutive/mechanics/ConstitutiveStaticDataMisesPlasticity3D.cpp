// $Id: ConstitutiveStaticDataBase.cpp 102 2009-11-11 10:47:23Z eckardt4 $

#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataMisesPlasticity3D.h"

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

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template<class Archive>
void NuTo::ConstitutiveStaticDataMisesPlasticity3D::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveStaticDataBase)
       & BOOST_SERIALIZATION_NVP(mEpsilonPEq)
       & BOOST_SERIALIZATION_NVP(mEpsilonP)
       & BOOST_SERIALIZATION_NVP(mSigmaB);
}
#endif // ENABLE_SERIALIZATION
