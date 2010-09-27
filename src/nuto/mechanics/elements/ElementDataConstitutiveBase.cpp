/*
 * ElementDataConstitutiveBase.cpp
 *
 *  Created on: Apr 16, 2010
 *      Author: unger3
 */

// $Id: $
#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/ConstitutiveBase.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/elements/ElementDataConstitutiveBase.h"
#include <assert.h>


NuTo::ElementDataConstitutiveBase::ElementDataConstitutiveBase() :  NuTo::ElementDataBase::ElementDataBase()
{
    //std::cout << "ElementDataConstitutiveBase constructor " << std::endl;
	mConstitutiveLaw = 0;
}
//! @brief deconstructor
NuTo::ElementDataConstitutiveBase::~ElementDataConstitutiveBase()
{
    //std::cout << "NuTo::ElementDataConstitutiveBase::~ElementDataConstitutiveBase()" << std::endl;
}

void NuTo::ElementDataConstitutiveBase::SetConstitutiveLaw(const ElementBase* rElement, NuTo::ConstitutiveBase* rConstitutiveLaw)
{
	//std::cout << "NuTo::ElementDataConstitutiveBase::SetConstitutiveLaw" << std::endl;
	assert(rConstitutiveLaw!=0);
	mConstitutiveLaw = rConstitutiveLaw;
	InitializeUpdatedConstitutiveLaw(rElement);
}

NuTo::ConstitutiveBase* NuTo::ElementDataConstitutiveBase::GetConstitutiveLaw(int rIp)
{
	return mConstitutiveLaw;
}

const NuTo::ConstitutiveBase* NuTo::ElementDataConstitutiveBase::GetConstitutiveLaw(int rIp)const
{
	return mConstitutiveLaw;
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::ElementDataConstitutiveBase::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ElementDataConstitutiveBase::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ElementDataConstitutiveBase::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ElementDataConstitutiveBase::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ElementDataConstitutiveBase::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ElementDataConstitutiveBase::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ElementDataConstitutiveBase::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ElementDataConstitutiveBase" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ElementDataBase)
       & BOOST_SERIALIZATION_NVP(mConstitutiveLaw);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ElementDataConstitutiveBase" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ElementDataConstitutiveBase)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::ElementDataConstitutiveBase)
#endif // ENABLE_SERIALIZATION
