// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/mechanics/GreenLagrangeStrain3D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient3D.h"
#include "nuto/mechanics/MechanicsException.h"

NuTo::GreenLagrangeStrain3D::GreenLagrangeStrain3D()
{
	mGreenLagrangeStrain[0] = 0.0;
	mGreenLagrangeStrain[1] = 0.0;
	mGreenLagrangeStrain[2] = 0.0;
	mGreenLagrangeStrain[3] = 0.0;
	mGreenLagrangeStrain[4] = 0.0;
	mGreenLagrangeStrain[5] = 0.0;
}

NuTo::GreenLagrangeStrain3D::GreenLagrangeStrain3D(const DeformationGradient3D& rDeformationGradient)
{
	throw MechanicsException("[NuTo::GreenLagrangeStrain3D::GreenLagrangeStrain3D] to be implemented.");
}

//! @brief ... get number of strain components
//! @return ... number of strain components
unsigned int NuTo::GreenLagrangeStrain3D::GetNumberOfComponents() const
{
    return 6;
}

//! @brief ... get Green Lagrange Strain
//! @return ... Green Strain (exx,eyy,ezz,gxy,gyz,gzx)
//! @sa mDeformationGradient
const double* NuTo::GreenLagrangeStrain3D::GetData() const
{
    return mGreenLagrangeStrain;
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::GreenLagrangeStrain3D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::GreenLagrangeStrain3D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::GreenLagrangeStrain3D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::GreenLagrangeStrain3D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::GreenLagrangeStrain3D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::GreenLagrangeStrain3D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::GreenLagrangeStrain3D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize GreenLagrangeStrain3D" << std::endl;
#endif
   ar & BOOST_SERIALIZATION_NVP(mGreenLagrangeStrain);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize GreenLagrangeStrain3D" << std::endl;
#endif
}
#endif // ENABLE_SERIALIZATION
