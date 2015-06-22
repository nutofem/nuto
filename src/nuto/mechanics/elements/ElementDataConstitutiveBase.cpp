// $Id$

/*
 * ElementDataConstitutiveBase.cpp
 *
 *  Created on: Apr 16, 2010
 *      Author: unger3
 */

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
	mConstitutiveLaw = nullptr;
}
//! @brief deconstructor
NuTo::ElementDataConstitutiveBase::~ElementDataConstitutiveBase()
{
    //std::cout << "NuTo::ElementDataConstitutiveBase::~ElementDataConstitutiveBase()" << std::endl;
}

void NuTo::ElementDataConstitutiveBase::SetConstitutiveLaw(const ElementBase* rElement, NuTo::ConstitutiveBase* rConstitutiveLaw)
{
	//std::cout << "NuTo::ElementDataConstitutiveBase::SetConstitutiveLaw" << std::endl;
	assert(rConstitutiveLaw!=nullptr);
    mConstitutiveLaw = rConstitutiveLaw;
	InitializeUpdatedConstitutiveLaw(rElement);
}

void NuTo::ElementDataConstitutiveBase::SetConstitutiveLaw(const ElementBase* rElement,int rIp, NuTo::ConstitutiveBase* rConstitutiveLaw)
{
	throw MechanicsException("[NuTo::ElementDataConstitutiveBase::SetConstitutiveLaw] one material for whole element.");
}

bool NuTo::ElementDataConstitutiveBase::HasConstitutiveLawAssigned(int rIp) const
{
	return mConstitutiveLaw != nullptr;
}


NuTo::ConstitutiveBase* NuTo::ElementDataConstitutiveBase::GetConstitutiveLaw(int rIp)
{
	if (mConstitutiveLaw==nullptr)
		throw MechanicsException("[NuTo::ElementDataConstitutiveBase::GetConstitutiveLaw] no constitutive law assigned yet.");
	return mConstitutiveLaw;
}

const NuTo::ConstitutiveBase* NuTo::ElementDataConstitutiveBase::GetConstitutiveLaw(int rIp)const
{
	if (mConstitutiveLaw==nullptr)
		throw MechanicsException("[NuTo::ElementDataConstitutiveBase::GetConstitutiveLaw] no constitutive law assigned yet.");
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
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::ElementDataConstitutiveBase)
#endif // ENABLE_SERIALIZATION
