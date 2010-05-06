/*
 * ElementDataConstitutiveBase.cpp
 *
 *  Created on: Apr 16, 2010
 *      Author: unger3
 */

// $Id: $

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/elements/ElementDataConstitutiveBase.h"
#include <assert.h>


NuTo::ElementDataConstitutiveBase::ElementDataConstitutiveBase() :  NuTo::ElementDataBase::ElementDataBase()
{
    std::cout << "ElementDataConstitutiveBase constructor " << std::endl;
	mConstitutiveLaw = 0;
}
//! @brief deconstructor
NuTo::ElementDataConstitutiveBase::~ElementDataConstitutiveBase()
{
    //std::cout << "NuTo::ElementDataConstitutiveBase::~ElementDataConstitutiveBase()" << std::endl;
}

void NuTo::ElementDataConstitutiveBase::SetConstitutiveLaw(const ElementWithDataBase* rElement, NuTo::ConstitutiveBase* rConstitutiveLaw)
{
	std::cout << "NuTo::ElementDataConstitutiveBase::SetConstitutiveLaw" << std::endl;
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
