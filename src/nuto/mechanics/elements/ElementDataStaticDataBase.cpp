// $Id: $
/*
 * ElementDataStaticDataBase.cpp
 *
 *  Created on: Apr 16, 2010
 *      Author: unger3
 */
#include "nuto/mechanics/elements/ElementDataStaticDataBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveBase.h"
#include <assert.h>


NuTo::ElementDataStaticDataBase::ElementDataStaticDataBase(const NuTo::IntegrationTypeBase* rIntegrationType) :
    NuTo::ElementDataBase::ElementDataBase(rIntegrationType)
{
}

NuTo::ElementDataStaticDataBase::~ElementDataStaticDataBase()
{
    for (int count=0; count< mStaticData.size(); count++)
    {
    	if (mStaticData[count]!=0)
    		delete mStaticData[count];
    }
}

NuTo::ConstitutiveStaticDataBase* NuTo::ElementDataStaticDataBase::GetStaticData(int rIp)
{
	assert(rIp<(int)mStaticData.size());
	return mStaticData[rIp];
}

const NuTo::ConstitutiveStaticDataBase* NuTo::ElementDataStaticDataBase::GetStaticData(int rIp)const
{
	assert(rIp<(int)mStaticData.size());
	return mStaticData[rIp];
}
