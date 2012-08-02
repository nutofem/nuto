// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain2D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient2D.h"

NuTo::EngineeringStrain2D::EngineeringStrain2D() : ConstitutiveOutputBase::ConstitutiveOutputBase()
{
	mEngineeringStrain[0] = 0.0;
	mEngineeringStrain[1] = 0.0;
	mEngineeringStrain[2] = 0.0;
}

NuTo::EngineeringStrain2D::EngineeringStrain2D(const DeformationGradient2D& rDeformationGradient)
{
    mEngineeringStrain[0] = rDeformationGradient.mDeformationGradient[0] -1;
    mEngineeringStrain[1] = rDeformationGradient.mDeformationGradient[3] -1;
    mEngineeringStrain[2] = rDeformationGradient.mDeformationGradient[1]+rDeformationGradient.mDeformationGradient[2];
}

NuTo::EngineeringStrain2D::EngineeringStrain2D(const EngineeringStrain2D& rEngineeringStrain)
{
    mEngineeringStrain[0] = rEngineeringStrain.mEngineeringStrain[0];
    mEngineeringStrain[1] = rEngineeringStrain.mEngineeringStrain[1];
    mEngineeringStrain[2] = rEngineeringStrain.mEngineeringStrain[2];
}

//! @brief ... get number of strain components
//! @return ... number of strain components
unsigned int NuTo::EngineeringStrain2D::GetNumberOfComponents() const
{
    return 3;
}

//! @brief ... get Engineering Strain
//! @return ... Engineering Strain (exx,eyy,gxy)
//! @sa mDeformationGradient
const double* NuTo::EngineeringStrain2D::GetData() const
{
    return mEngineeringStrain;
}

//! @brief ... set Engineering Strain
//! @return ... Engineering Strain (exx,eyy,gxy)
//! @sa mDeformationGradient
void NuTo::EngineeringStrain2D::SetData(double rEngineeringStrain[3])
{
	mEngineeringStrain[0] = rEngineeringStrain[0];
	mEngineeringStrain[1] = rEngineeringStrain[1];
	mEngineeringStrain[2] = rEngineeringStrain[2];
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::EngineeringStrain2D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::EngineeringStrain2D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::EngineeringStrain2D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::EngineeringStrain2D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::EngineeringStrain2D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::EngineeringStrain2D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::EngineeringStrain2D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize EngineeringStrain2D" << std::endl;
#endif
    ar &  BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveOutputBase)
       & BOOST_SERIALIZATION_NVP(mEngineeringStrain);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize EngineeringStrain2D" << std::endl;
#endif
}
#endif // ENABLE_SERIALIZATION
