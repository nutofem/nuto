// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/mechanics/LocalEqPlasticStrain.h"

NuTo::LocalEqPlasticStrain::LocalEqPlasticStrain(): ConstitutiveOutputBase::ConstitutiveOutputBase(), FullVector<double,2>()
{
}


#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::LocalEqPlasticStrain::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::LocalEqPlasticStrain::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::LocalEqPlasticStrain::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::LocalEqPlasticStrain::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::LocalEqPlasticStrain::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::LocalEqPlasticStrain::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::LocalEqPlasticStrain::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize LocalEqPlasticStrain" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveOutputBase);
    ar & boost::serialization::make_nvp ("LocalEqPlasticStrainEigen",boost::serialization::base_object< FullVector<double,2> > ( *this ));
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize LocalEqPlasticStrain" << std::endl;
#endif
}
#endif // ENABLE_SERIALIZATION

