// $Id: ConstitutiveStaticDataBase.cpp 102 2009-11-11 10:47:23Z eckardt4 $

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

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
template void NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain2DPlaneStrain::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain2DPlaneStrain::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain2DPlaneStrain::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain2DPlaneStrain::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain2DPlaneStrain::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain2DPlaneStrain::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain2DPlaneStrain::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstitutiveStaticDataPrevEngineeringStressStrain2DPlaneStrain" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveStaticDataBase)
       & BOOST_SERIALIZATION_NVP(mPrevSigma)
       & BOOST_SERIALIZATION_NVP(mPrevStrain)
       & BOOST_SERIALIZATION_NVP(mPrevTotalEnergy)
       & BOOST_SERIALIZATION_NVP(mPrevElasticEnergy);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstitutiveStaticDataPrevEngineeringStressStrain2DPlaneStrain" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain2DPlaneStrain)
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


