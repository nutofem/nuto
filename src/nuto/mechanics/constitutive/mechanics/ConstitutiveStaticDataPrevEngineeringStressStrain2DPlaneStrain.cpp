// $Id: ConstitutiveStaticDataBase.cpp 102 2009-11-11 10:47:23Z eckardt4 $

#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataPrevEngineeringStressStrain2DPlaneStrain.h"

//! @brief constructor
NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain2DPlaneStrain::ConstitutiveStaticDataPrevEngineeringStressStrain2DPlaneStrain() : ConstitutiveStaticDataBase()
{
    mPrevSigma[0] = 0.;
    mPrevSigma[1] = 0.;
    mPrevSigma[2] = 0.;
    mPrevSigma[3] = 0.;

    mPrevStrain[0] = 0.;
    mPrevStrain[1] = 0.;
    mPrevStrain[2] = 0.;

    mPrevTotalEnergy = 0.;

    mPrevElasticEnergy = 0.;
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template<class Archive>
void NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain2DPlaneStrain::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveStaticDataBase)
       & BOOST_SERIALIZATION_NVP(mPrevSigma)
       & BOOST_SERIALIZATION_NVP(mPrevStrain)
       & BOOST_SERIALIZATION_NVP(mPrevTotalEnergy)
       & BOOST_SERIALIZATION_NVP(mPrevElasticEnergy);
}
#endif // ENABLE_SERIALIZATION

//! brief set the previous stress
void NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain2DPlaneStrain::SetPrevStress(const double rPrevSigma[4])
{
	mPrevSigma[0] = rPrevSigma[0];
	mPrevSigma[1] = rPrevSigma[1];
	mPrevSigma[2] = rPrevSigma[2];
	mPrevSigma[3] = rPrevSigma[3];
}

//! brief set the previous stress
void NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain2DPlaneStrain::SetPrevStrain(const double rPrevStrain[3])
{
	mPrevStrain[0] = rPrevStrain[0];
	mPrevStrain[1] = rPrevStrain[1];
	mPrevStrain[2] = rPrevStrain[2];
}


