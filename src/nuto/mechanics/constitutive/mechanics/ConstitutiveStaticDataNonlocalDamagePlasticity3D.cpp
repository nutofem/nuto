// $ld: $ 
// ConstitutiveStaticDataNonlocalDamagePlasticity3D.cpp
// created May 6, 2010 by Joerg F. Unger


#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataNonlocalDamagePlasticity3D.h"

//! @brief constructor
NuTo::ConstitutiveStaticDataNonlocalDamagePlasticity3D::ConstitutiveStaticDataNonlocalDamagePlasticity3D() : ConstitutiveStaticDataPrevEngineeringStressStrain3D()
{
	mOmega = 0.;
	mKappa = 0.;

	mEpsilonP[0] = 0.;
	mEpsilonP[1] = 0.;
	mEpsilonP[2] = 0.;
	mEpsilonP[3] = 0.;
	mEpsilonP[4] = 0.;
	mEpsilonP[5] = 0.;
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template<class Archive>
void NuTo::ConstitutiveStaticDataNonlocalDamagePlasticity3D::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveStaticDataPrevEngineeringStressStrain3D)
       & BOOST_SERIALIZATION_NVP(mOmega)
       & BOOST_SERIALIZATION_NVP(mKappa)
       & BOOST_SERIALIZATION_NVP(mEpsilonP);
}
#endif // ENABLE_SERIALIZATION
