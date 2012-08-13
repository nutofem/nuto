// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain3D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient3D.h"

NuTo::EngineeringStrain3D::EngineeringStrain3D() : ConstitutiveOutputBase::ConstitutiveOutputBase()
{
	mEngineeringStrain[0] = 0.0;
	mEngineeringStrain[1] = 0.0;
	mEngineeringStrain[2] = 0.0;
	mEngineeringStrain[3] = 0.0;
	mEngineeringStrain[4] = 0.0;
	mEngineeringStrain[5] = 0.0;
}

NuTo::EngineeringStrain3D::EngineeringStrain3D(const DeformationGradient3D& rDeformationGradient)
{
    mEngineeringStrain[0] = rDeformationGradient.mDeformationGradient[0] -1.;
    mEngineeringStrain[1] = rDeformationGradient.mDeformationGradient[4] -1.;
    mEngineeringStrain[2] = rDeformationGradient.mDeformationGradient[8] -1.;
    mEngineeringStrain[3] = rDeformationGradient.mDeformationGradient[1]+rDeformationGradient.mDeformationGradient[3];
    mEngineeringStrain[4] = rDeformationGradient.mDeformationGradient[5]+rDeformationGradient.mDeformationGradient[7];
    mEngineeringStrain[5] = rDeformationGradient.mDeformationGradient[2]+rDeformationGradient.mDeformationGradient[6];
}

//! @brief ... get number of strain components
//! @return ... number of strain components
unsigned int NuTo::EngineeringStrain3D::GetNumberOfComponents() const
{
    return 6;
}

//! @brief ... get Engineering Strain
//! @return ... Engineering Strain (exx,eyy,ezz,gxy,gyz,gzx)
//! @sa mDeformationGradient
const double* NuTo::EngineeringStrain3D::GetData() const
{
    return mEngineeringStrain;
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::EngineeringStrain3D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::EngineeringStrain3D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::EngineeringStrain3D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::EngineeringStrain3D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::EngineeringStrain3D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::EngineeringStrain3D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::EngineeringStrain3D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize EngineeringStrain3D" << std::endl;
#endif
    ar &  BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveOutputBase)
       & BOOST_SERIALIZATION_NVP(mEngineeringStrain);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize EngineeringStrain3D" << std::endl;
#endif
}
#endif // ENABLE_SERIALIZATION
