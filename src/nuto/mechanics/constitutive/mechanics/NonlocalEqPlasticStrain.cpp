// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/mechanics/NonlocalEqPlasticStrain.h"

NuTo::NonlocalEqPlasticStrain::NonlocalEqPlasticStrain() : ConstitutiveInputBase::ConstitutiveInputBase(), FullVector<double,2>()
{
}


#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::NonlocalEqPlasticStrain::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::NonlocalEqPlasticStrain::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::NonlocalEqPlasticStrain::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::NonlocalEqPlasticStrain::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::NonlocalEqPlasticStrain::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::NonlocalEqPlasticStrain::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::NonlocalEqPlasticStrain::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize NonlocalEqPlasticStrain" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveInputBase);
    ar & boost::serialization::make_nvp ("NonlocalEqPlasticStrainEigen",boost::serialization::base_object< FullVector<double,2> > ( *this ));
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize NonlocalEqPlasticStrain" << std::endl;
#endif
}
#endif // ENABLE_SERIALIZATION

