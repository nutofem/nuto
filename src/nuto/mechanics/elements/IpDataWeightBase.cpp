// $Id: IpDataWeightBase.cpp 569 2011-08-16 21:13:55Z unger3 $
// IpDataWeightBase.cpp
// created Apr 29, 2010 by Joerg F. Unger

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/elements/IpDataWeightBase.h"
#include "nuto/mechanics/MechanicsException.h"

NuTo::IpDataWeightBase::IpDataWeightBase() : IpDataBase()
{
	mWeight=0.;
}

NuTo::IpDataWeightBase::~IpDataWeightBase()
{
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::IpDataWeightBase::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::IpDataWeightBase::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::IpDataWeightBase::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::IpDataWeightBase::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::IpDataWeightBase::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::IpDataWeightBase::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::IpDataWeightBase::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize IpDataWeightBase" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(IpDataBase)
       & BOOST_SERIALIZATION_NVP(mWeight);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize IpDataWeightBase" << std::endl;
#endif
}
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::IpDataWeightBase)
#endif // ENABLE_SERIALIZATION
