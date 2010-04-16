// $Id: ConstitutiveStaticDataBase.cpp 102 2009-11-11 10:47:23Z eckardt4 $

#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataPrevEngineeringStressStrain3D.h"

//! @brief constructor
NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain3D::ConstitutiveStaticDataPrevEngineeringStressStrain3D() : ConstitutiveStaticDataBase()
{
    mPrevSigma[0] = 0.;
    mPrevSigma[1] = 0.;
    mPrevSigma[2] = 0.;
    mPrevSigma[3] = 0.;
    mPrevSigma[4] = 0.;
    mPrevSigma[5] = 0.;

    mPrevStrain[0] = 0.;
    mPrevStrain[1] = 0.;
    mPrevStrain[2] = 0.;
    mPrevStrain[3] = 0.;
    mPrevStrain[4] = 0.;
    mPrevStrain[5] = 0.;

    mPrevTotalEnergy = 0.;

    mPrevElasticEnergy = 0.;
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template<class Archive>
void NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain3D::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveStaticDataBase)
       & BOOST_SERIALIZATION_NVP(mPrevSigma)
       & BOOST_SERIALIZATION_NVP(mPrevStrain)
       & BOOST_SERIALIZATION_NVP(mPrevTotalEnergy)
       & BOOST_SERIALIZATION_NVP(mPrevElasticEnergy);
}
#endif // ENABLE_SERIALIZATION

//! brief set the previous stress
void NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain3D::SetPrevStress(const double rPrevSigma[6])
{
	mPrevSigma[0] = rPrevSigma[0];
	mPrevSigma[1] = rPrevSigma[1];
	mPrevSigma[2] = rPrevSigma[2];
	mPrevSigma[3] = rPrevSigma[3];
	mPrevSigma[4] = rPrevSigma[4];
	mPrevSigma[5] = rPrevSigma[5];
}

//! brief set the previous stress
void NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain3D::SetPrevStrain(const double rPrevStrain[6])
{
	mPrevStrain[0] = rPrevStrain[0];
	mPrevStrain[1] = rPrevStrain[1];
	mPrevStrain[2] = rPrevStrain[2];
	mPrevStrain[3] = rPrevStrain[3];
	mPrevStrain[4] = rPrevStrain[4];
	mPrevStrain[5] = rPrevStrain[5];
}


