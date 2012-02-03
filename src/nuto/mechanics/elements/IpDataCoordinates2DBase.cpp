// $Id: IpDataCoordinates2DBase.cpp 569 2011-08-16 21:13:55Z unger3 $
// IpDataCoordinates2DBase.cpp
// created Apr 29, 2010 by Joerg F. Unger

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/elements/IpDataCoordinates2DBase.h"
#include "nuto/mechanics/MechanicsException.h"

NuTo::IpDataCoordinates2DBase::IpDataCoordinates2DBase() : IpDataBase()
{
	mLocalCoordinates[0]=0.;
	mLocalCoordinates[1]=0.;
}

NuTo::IpDataCoordinates2DBase::~IpDataCoordinates2DBase()
{
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::IpDataCoordinates2DBase::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::IpDataCoordinates2DBase::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::IpDataCoordinates2DBase::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::IpDataCoordinates2DBase::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::IpDataCoordinates2DBase::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::IpDataCoordinates2DBase::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::IpDataCoordinates2DBase::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize IpDataCoordinates2DBase" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(IpDataBase)
       & BOOST_SERIALIZATION_NVP(mLocalCoordinates);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize IpDataCoordinates2DBase" << std::endl;
#endif
}
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::IpDataCoordinates2DBase)
#endif // ENABLE_SERIALIZATION
