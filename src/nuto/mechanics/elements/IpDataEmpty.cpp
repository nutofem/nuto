// $Id$ 
// IpDataNonlocalBase.cpp
// created Apr 28, 2010 by Joerg F. Unger

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include <assert.h>
#include <iostream>
#include "nuto/mechanics/elements/IpDataEmpty.h"

NuTo::IpDataEmpty::IpDataEmpty()
{
	//std::cout << std::endl << "call Consctructor [NuTo::IpDataEmpty::~IpDataEmpty]." << std::endl;
}

void NuTo::IpDataEmpty::Initialize(const ElementBase* rElement, const ConstitutiveBase* rConstitutive)
{
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::IpDataEmpty::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::IpDataEmpty::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::IpDataEmpty::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::IpDataEmpty::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::IpDataEmpty::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::IpDataEmpty::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::IpDataEmpty::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize IpDataEmpty" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(IpDataBase);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize IpDataEmpty" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::IpDataEmpty)
#endif // ENABLE_SERIALIZATION
