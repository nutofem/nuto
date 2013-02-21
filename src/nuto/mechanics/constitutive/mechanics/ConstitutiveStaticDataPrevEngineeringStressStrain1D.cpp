// $Id: ConstitutiveStaticDataBase.cpp 102 2009-11-11 10:47:23Z eckardt4 $

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataPrevEngineeringStressStrain1D.h"

//! @brief constructor
NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain1D::ConstitutiveStaticDataPrevEngineeringStressStrain1D() : ConstitutiveStaticDataBase()
{
    mPrevTotalEnergy = 0.;
    mPrevElasticEnergy = 0.;
}

NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain1D& NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain1D::operator= (NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain1D const& rOther)
{
    NuTo::ConstitutiveStaticDataBase::operator= (rOther);
    mPrevSigma = rOther.mPrevSigma;
    mPrevStrain = rOther.mPrevStrain;
    mPrevTotalEnergy = rOther.mPrevTotalEnergy;
    mPrevElasticEnergy = rOther.mPrevElasticEnergy;
    return *this;
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain1D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain1D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain1D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain1D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain1D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain1D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain1D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstitutiveStaticDataPrevEngineeringStressStrain1D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveStaticDataBase)
       & BOOST_SERIALIZATION_NVP(mPrevSigma)
       & BOOST_SERIALIZATION_NVP(mPrevStrain)
       & BOOST_SERIALIZATION_NVP(mPrevTotalEnergy)
       & BOOST_SERIALIZATION_NVP(mPrevElasticEnergy);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstitutiveStaticDataPrevEngineeringStressStrain1D" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain1D)
#endif // ENABLE_SERIALIZATION

//! brief set the previous stress
void NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain1D::SetPrevStress(const EngineeringStress1D& rPrevSigma)
{
	mPrevSigma = rPrevSigma;
}

//! brief set the previous stress
void NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain1D::SetPrevStrain(const EngineeringStrain1D& rPrevStrain)
{
	mPrevStrain = rPrevStrain;
}


