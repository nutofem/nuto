// $Id: IpDataStaticDataWeightCoordinates2D.cpp 569 2011-08-16 21:13:55Z unger3 $
// IpDataStaticDataWeightCoordinates2D.cpp
// created Apr 29, 2010 by Joerg F. Unger

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/ConstitutiveStaticDataBase.h"
#include "nuto/mechanics/elements/IpDataStaticDataWeightCoordinates2D.h"
#include "nuto/mechanics/elements/ElementBase.h"
NuTo::IpDataStaticDataWeightCoordinates2D::IpDataStaticDataWeightCoordinates2D() :NuTo::IpDataBase::IpDataBase() ,
    NuTo::IpDataStaticDataBase::IpDataStaticDataBase() ,
    NuTo::IpDataCoordinates2DBase::IpDataCoordinates2DBase(),
    NuTo::IpDataWeightBase::IpDataWeightBase()
{

}

NuTo::IpDataStaticDataWeightCoordinates2D::~IpDataStaticDataWeightCoordinates2D()
{
}

void NuTo::IpDataStaticDataWeightCoordinates2D::Initialize(const ElementBase* rElement, const ConstitutiveBase* rConstitutive)
{
	if (mStaticData!=0)
		delete mStaticData;
	if (rConstitutive!=0)
	    mStaticData = rElement->AllocateStaticData(rConstitutive);
	else
		mStaticData = 0;
}

//! @brief returns the enum of IP data type
//! @return enum of IPDataType
NuTo::IpData::eIpDataType NuTo::IpDataStaticDataWeightCoordinates2D::GetIpDataType()const
{
    return NuTo::IpData::STATICDATAWEIGHTCOORDINATES2D;
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::IpDataStaticDataWeightCoordinates2D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::IpDataStaticDataWeightCoordinates2D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::IpDataStaticDataWeightCoordinates2D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::IpDataStaticDataWeightCoordinates2D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::IpDataStaticDataWeightCoordinates2D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::IpDataStaticDataWeightCoordinates2D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::IpDataStaticDataWeightCoordinates2D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize IpDataStaticDataWeightCoordinates2D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(IpDataStaticDataBase)
       & BOOST_SERIALIZATION_BASE_OBJECT_NVP(IpDataCoordinates2DBase)
    & BOOST_SERIALIZATION_BASE_OBJECT_NVP(IpDataWeightBase);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize IpDataStaticDataWeightCoordinates2D" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::IpDataStaticDataWeightCoordinates2D)
#endif // ENABLE_SERIALIZATION
