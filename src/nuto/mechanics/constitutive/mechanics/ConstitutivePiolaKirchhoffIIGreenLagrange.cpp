// $Id: ConstitutivePiolaKirchhoffIIGreenLagrange.cpp 112 2009-11-17 16:43:15Z unger3 $

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/constitutive/mechanics/ConstitutivePiolaKirchhoffIIGreenLagrange.h"

NuTo::ConstitutivePiolaKirchhoffIIGreenLagrange::ConstitutivePiolaKirchhoffIIGreenLagrange() : ConstitutiveBase()
{
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::ConstitutivePiolaKirchhoffIIGreenLagrange::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutivePiolaKirchhoffIIGreenLagrange::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutivePiolaKirchhoffIIGreenLagrange::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutivePiolaKirchhoffIIGreenLagrange::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutivePiolaKirchhoffIIGreenLagrange::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutivePiolaKirchhoffIIGreenLagrange::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstitutivePiolaKirchhoffIIGreenLagrange::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstitutiveEngineeringStressStrain" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveBase);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstitutiveEngineeringStressStrain" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstitutivePiolaKirchhoffIIGreenLagrange)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::ConstitutivePiolaKirchhoffIIGreenLagrange)
#endif // ENABLE_SERIALIZATION

