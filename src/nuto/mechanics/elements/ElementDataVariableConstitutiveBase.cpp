// $Id$

/*
 * ElementDataVariableConstitutiveBase.cpp
 *
 *  Created on: Oct 30, 2014
 *      Author: keszler
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
#include "nuto/mechanics/elements/ElementDataVariableConstitutiveBase.h"
#include <assert.h>


NuTo::ElementDataVariableConstitutiveBase::ElementDataVariableConstitutiveBase() :  NuTo::ElementDataBase::ElementDataBase()
{
    //std::cout << "ElementDataVariableConstitutiveBase constructor " << std::endl;
	//mVarConstitutiveLaw.clear();
}
//! @brief deconstructor
NuTo::ElementDataVariableConstitutiveBase::~ElementDataVariableConstitutiveBase()
{
//	boost::ptr_vector<ConstitutiveBase>::iterator it=mVarConstitutiveLaw.begin();
//	for(size_t i=0;i< mVarConstitutiveLaw.size();++i)
//	{
//		mVarConstitutiveLaw.insert(it,0);
//		++it;
//	}
//	std::cout << "NuTo::ElementDataVariableConstitutiveBase::~ElementDataVariableConstitutiveBase()" << std::endl;
}

void NuTo::ElementDataVariableConstitutiveBase::SetConstitutiveLaw(const ElementBase* rElement, NuTo::ConstitutiveBase* rConstitutiveLaw)
{
	throw MechanicsException("[NuTo::ElementDataVariableConstitutiveBase::SetConstitutiveLaw] needs different materials for one element.");
}


void NuTo::ElementDataVariableConstitutiveBase::SetConstitutiveLaw(const ElementBase* rElement, int rIp, NuTo::ConstitutiveBase* rConstitutiveLaw)
{
	//std::cout << "NuTo::ElementDataVariableConstitutiveBase::SetConstitutiveLaw" << std::endl;
	assert(rConstitutiveLaw!=0);

	if (rIp>=(int) mVarConstitutiveLaw.size())
		mVarConstitutiveLaw.push_back(rConstitutiveLaw);
	else
	{
		mVarConstitutiveLaw[rIp]=rConstitutiveLaw;

	}

	InitializeUpdatedConstitutiveLaw(rElement,rIp);
}

bool NuTo::ElementDataVariableConstitutiveBase::HasConstitutiveLawAssigned(int rIp)
{
	bool check=true;
	if(mVarConstitutiveLaw[rIp]==0)
		check=!check;
	return check;
}


NuTo::ConstitutiveBase* NuTo::ElementDataVariableConstitutiveBase::GetConstitutiveLaw(int rIp)
{
	if(mVarConstitutiveLaw[rIp]==0)
		throw MechanicsException("[NuTo::ElementDataVariableConstitutiveBase::GetConstitutiveLaw] no constitutive law assigned yet.");
	return mVarConstitutiveLaw[rIp];
	//return mVarConstitutiveLaw[rIp];
}

const NuTo::ConstitutiveBase* NuTo::ElementDataVariableConstitutiveBase::GetConstitutiveLaw(int rIp)const
{
	if(mVarConstitutiveLaw[rIp]==0)
		throw MechanicsException("[NuTo::ElementDataVariableConstitutiveBase::GetConstitutiveLaw] no constitutive law assigned yet.");
	return mVarConstitutiveLaw[rIp];
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::ElementDataVariableConstitutiveBase::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ElementDataVariableConstitutiveBase::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ElementDataVariableConstitutiveBase::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ElementDataVariableConstitutiveBase::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ElementDataVariableConstitutiveBase::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ElementDataVariableConstitutiveBase::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ElementDataVariableConstitutiveBase::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ElementDataVariableConstitutiveBase" << std::endl;
#endif
//    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ElementDataBase)
//       & BOOST_SERIALIZATION_NVP(mVarConstitutiveLaw);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ElementDataVariableConstitutiveBase" << std::endl;
#endif
}
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::ElementDataVariableConstitutiveBase)
#endif // ENABLE_SERIALIZATION
