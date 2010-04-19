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


NuTo::ElementDataConstitutiveBase::ElementDataConstitutiveBase(const NuTo::IntegrationTypeBase* rIntegrationType) :  NuTo::ElementDataBase::ElementDataBase(rIntegrationType)
{
    //std::cout << "ElementDataConstitutiveBase constructor " << std::endl;
	mConstitutiveLaw = 0;
}

void NuTo::ElementDataConstitutiveBase::SetConstitutiveLaw(const ElementBase* rElement, NuTo::ConstitutiveBase* rConstitutiveLaw)
{
	assert(rConstitutiveLaw!=0);
	mConstitutiveLaw = rConstitutiveLaw;
	UpdateForModifiedConstitutiveLaw(rElement);
}

NuTo::ConstitutiveBase* NuTo::ElementDataConstitutiveBase::GetConstitutiveLaw(int rIp)
{
	assert(mConstitutiveLaw!=0);
	return mConstitutiveLaw;
}

const NuTo::ConstitutiveBase* NuTo::ElementDataConstitutiveBase::GetConstitutiveLaw(int rIp)const
{
	assert(mConstitutiveLaw!=0);
	return mConstitutiveLaw;
}
