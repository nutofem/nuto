// $Id: ConstitutiveStaticDataBase.cpp 102 2009-11-11 10:47:23Z eckardt4 $

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

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

NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain3D& NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain3D::operator= (NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain3D const& rOther)
{
    NuTo::ConstitutiveStaticDataBase::operator= (rOther);
    mPrevSigma[0] = rOther.mPrevSigma[0];
    mPrevSigma[1] = rOther.mPrevSigma[1];
    mPrevSigma[2] = rOther.mPrevSigma[2];
    mPrevSigma[3] = rOther.mPrevSigma[3];
    mPrevSigma[4] = rOther.mPrevSigma[4];
    mPrevSigma[5] = rOther.mPrevSigma[5];
    mPrevStrain[0] = rOther.mPrevStrain[0];
    mPrevStrain[1] = rOther.mPrevStrain[1];
    mPrevStrain[2] = rOther.mPrevStrain[2];
    mPrevStrain[3] = rOther.mPrevStrain[3];
    mPrevStrain[4] = rOther.mPrevStrain[4];
    mPrevStrain[5] = rOther.mPrevStrain[5];
    mPrevTotalEnergy = rOther.mPrevTotalEnergy;
    mPrevElasticEnergy = rOther.mPrevElasticEnergy;
    return *this;
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain3D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain3D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain3D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain3D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain3D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain3D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain3D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstitutiveStaticDataPrevEngineeringStressStrain3D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveStaticDataBase)
       & BOOST_SERIALIZATION_NVP(mPrevSigma)
       & BOOST_SERIALIZATION_NVP(mPrevStrain)
       & BOOST_SERIALIZATION_NVP(mPrevTotalEnergy)
       & BOOST_SERIALIZATION_NVP(mPrevElasticEnergy);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstitutiveStaticDataPrevEngineeringStressStrain3D" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain3D)
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


