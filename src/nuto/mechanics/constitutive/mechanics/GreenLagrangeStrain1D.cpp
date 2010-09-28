
#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/mechanics/GreenLagrangeStrain1D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient1D.h"

NuTo::GreenLagrangeStrain1D::GreenLagrangeStrain1D()
{
	mGreenLagrangeStrain = 0.0;
}

NuTo::GreenLagrangeStrain1D::GreenLagrangeStrain1D(const DeformationGradient1D& rDeformationGradient)
{
    mGreenLagrangeStrain = 0.5*(rDeformationGradient.GetDeformationGradient1D()[0]*rDeformationGradient.GetDeformationGradient1D()[0] -1);
}

//! @brief ... get number of strain components
//! @return ... number of strain components
unsigned int NuTo::GreenLagrangeStrain1D::GetNumberOfComponents() const
{
    return 1;
}

//! @brief ... get Green Strain
//! @return ... Green Strain (exx)
const double* NuTo::GreenLagrangeStrain1D::GetData() const
{
    return &mGreenLagrangeStrain;
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::GreenLagrangeStrain1D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::GreenLagrangeStrain1D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::GreenLagrangeStrain1D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::GreenLagrangeStrain1D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::GreenLagrangeStrain1D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::GreenLagrangeStrain1D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::GreenLagrangeStrain1D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize GreenLagrangeStrain1D" << std::endl;
#endif
   ar & BOOST_SERIALIZATION_NVP(mGreenLagrangeStrain);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize GreenLagrangeStrain1D" << std::endl;
#endif
}
#endif // ENABLE_SERIALIZATION
