// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/thermal/TemperatureGradient3D.h"

NuTo::TemperatureGradient3D::TemperatureGradient3D() : ConstitutiveInputBase()
{
	mTemperatureGradient[0] = 0.0;
	mTemperatureGradient[1] = 0.0;
	mTemperatureGradient[2] = 0.0;
}

//! @brief ... get number of strain components
//! @return ... number of strain components
unsigned int NuTo::TemperatureGradient3D::GetNumberOfComponents() const
{
    return 3;
}

//! @brief ... get Engineering Strain
//! @return ... Engineering Strain (exx)
const double* NuTo::TemperatureGradient3D::GetData() const
{
    return mTemperatureGradient;
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::TemperatureGradient3D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::TemperatureGradient3D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::TemperatureGradient3D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::TemperatureGradient3D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::TemperatureGradient3D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::TemperatureGradient3D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::TemperatureGradient3D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize TemperatureGradient3D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveInputBase)
       & BOOST_SERIALIZATION_NVP(mTemperatureGradient);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize TemperatureGradient3D" << std::endl;
#endif
}
#endif // ENABLE_SERIALIZATION

