
#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain1D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient1D.h"

NuTo::EngineeringStrain1D::EngineeringStrain1D()
{
	mEngineeringStrain = 0.0;
}

NuTo::EngineeringStrain1D::EngineeringStrain1D(const DeformationGradient1D& rDeformationGradient)
{
    mEngineeringStrain = rDeformationGradient.GetDeformationGradient1D()[0] -1;
}

//! @brief ... get number of strain components
//! @return ... number of strain components
unsigned int NuTo::EngineeringStrain1D::GetNumberOfComponents() const
{
    return 1;
}

//! @brief ... get Engineering Strain
//! @return ... Engineering Strain (exx)
const double* NuTo::EngineeringStrain1D::GetData() const
{
    return &mEngineeringStrain;
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::EngineeringStrain1D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::EngineeringStrain1D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::EngineeringStrain1D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::EngineeringStrain1D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::EngineeringStrain1D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::EngineeringStrain1D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::EngineeringStrain1D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize EngineeringStrain1D" << std::endl;
#endif
   ar & BOOST_SERIALIZATION_NVP(mEngineeringStrain);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize EngineeringStrain1D" << std::endl;
#endif
}
#endif // ENABLE_SERIALIZATION

