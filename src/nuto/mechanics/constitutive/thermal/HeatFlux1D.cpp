// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/thermal/HeatFlux1D.h"

NuTo::HeatFlux1D::HeatFlux1D(): ConstitutiveOutputBase()
{
	mHeatFlux = 0.0;
}


//! @brief ... get number of strain components
//! @return ... number of strain components
unsigned int NuTo::HeatFlux1D::GetNumberOfComponents() const
{
    return 1;
}

//! @brief ... get Engineering Strain
//! @return ... Engineering Strain (exx)
const double* NuTo::HeatFlux1D::GetData() const
{
    return &mHeatFlux;
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::HeatFlux1D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::HeatFlux1D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::HeatFlux1D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::HeatFlux1D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::HeatFlux1D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::HeatFlux1D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::HeatFlux1D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize HeatFlux1D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveOutputBase)
       & BOOST_SERIALIZATION_NVP(mHeatFlux);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize HeatFlux1D" << std::endl;
#endif
}
#endif // ENABLE_SERIALIZATION

