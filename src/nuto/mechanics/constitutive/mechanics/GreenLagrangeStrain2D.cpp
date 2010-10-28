// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/mechanics/GreenLagrangeStrain2D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient2D.h"

NuTo::GreenLagrangeStrain2D::GreenLagrangeStrain2D()
{
	mGreenLagrangeStrain[0] = 0.0;
	mGreenLagrangeStrain[1] = 0.0;
	mGreenLagrangeStrain[2] = 0.0;
}


NuTo::GreenLagrangeStrain2D::GreenLagrangeStrain2D(const DeformationGradient2D& rDeformationGradient)
{
    mGreenLagrangeStrain[0] = 0.5 * (rDeformationGradient.GetDeformationGradient2D()[0] * rDeformationGradient.GetDeformationGradient2D()[0] + rDeformationGradient.GetDeformationGradient2D()[1] * rDeformationGradient.GetDeformationGradient2D()[1] - 1.0);
    mGreenLagrangeStrain[1] = 0.5 * (rDeformationGradient.GetDeformationGradient2D()[3] * rDeformationGradient.GetDeformationGradient2D()[3] + rDeformationGradient.GetDeformationGradient2D()[2] * rDeformationGradient.GetDeformationGradient2D()[2] - 1.0);
    mGreenLagrangeStrain[2] = rDeformationGradient.GetDeformationGradient2D()[0] * rDeformationGradient.GetDeformationGradient2D()[2] + rDeformationGradient.GetDeformationGradient2D()[1] * rDeformationGradient.GetDeformationGradient2D()[3];
}

//! @brief ... get number of strain components
//! @return ... number of strain components
unsigned int NuTo::GreenLagrangeStrain2D::GetNumberOfComponents() const
{
    return 3;
}

//! @brief ... get green lagrange Strain
//! @return ... (exx,eyy,gxy)
//! @sa mDeformationGradient
const double* NuTo::GreenLagrangeStrain2D::GetData() const
{
    return mGreenLagrangeStrain;
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::GreenLagrangeStrain2D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::GreenLagrangeStrain2D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::GreenLagrangeStrain2D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::GreenLagrangeStrain2D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::GreenLagrangeStrain2D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::GreenLagrangeStrain2D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::GreenLagrangeStrain2D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize GreenLagrangeStrain2D" << std::endl;
#endif
   ar & BOOST_SERIALIZATION_NVP(mGreenLagrangeStrain);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize GreenLagrangeStrain2D" << std::endl;
#endif
}
#endif // ENABLE_SERIALIZATION
