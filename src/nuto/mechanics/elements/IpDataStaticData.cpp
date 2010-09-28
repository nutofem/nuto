// $ld: $ 
// IpDataStaticData.cpp
// created Apr 29, 2010 by Joerg F. Unger

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/elements/IpDataStaticData.h"
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveStaticDataBase.h"

NuTo::IpDataStaticData::IpDataStaticData() : IpDataBase()
{
	mStaticData = 0;
}

NuTo::IpDataStaticData::~IpDataStaticData()
{
}

void NuTo::IpDataStaticData::Initialize(const ElementBase* rElement, const ConstitutiveBase* rConstitutive)
{
	if (mStaticData!=0)
		delete mStaticData;
	if (rConstitutive!=0)
	    mStaticData = rElement->AllocateStaticData(rConstitutive);
	else
		mStaticData = 0;
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::IpDataStaticData::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::IpDataStaticData::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::IpDataStaticData::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::IpDataStaticData::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::IpDataStaticData::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::IpDataStaticData::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::IpDataStaticData::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize IpDataStaticData" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(IpDataStaticDataBase);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize IpDataStaticData" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::IpDataStaticData)
#endif // ENABLE_SERIALIZATION
