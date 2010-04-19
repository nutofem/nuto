// $Id: $
#include "nuto/mechanics/elements/ElementDataConstitutiveStaticData.h"
#include "nuto/mechanics/constitutive/ConstitutiveBase.h"
#include <assert.h>


NuTo::ElementDataConstitutiveStaticData::ElementDataConstitutiveStaticData(const NuTo::IntegrationTypeBase* rIntegrationType) :
ElementDataConstitutiveBase(rIntegrationType) , ElementDataStaticDataBase(rIntegrationType), NuTo::ElementDataBase::ElementDataBase(rIntegrationType)
{
}

NuTo::ElementDataConstitutiveStaticData::~ElementDataConstitutiveStaticData()
{
}

//! @brief update the information related to a modification of the integration type, e.g. reallocation of the static data
void NuTo::ElementDataConstitutiveStaticData::UpdateForModifiedIntegrationType(const ElementBase* rElement)
{
	mStaticData.resize(mIntegrationType->GetNumIntegrationPoints());
	if (mConstitutiveLaw!=0)
	{
		for (unsigned int count=0; count<mStaticData.size(); count++)
		{
			mStaticData[count] = rElement->AllocateStaticData(mConstitutiveLaw);
		}
	}
}

//! @brief update the information related to a modification of the constitutive law, e.g. reallocation of the static data
void NuTo::ElementDataConstitutiveStaticData::UpdateForModifiedConstitutiveLaw(const ElementBase* rElement)
{
	mStaticData.resize(mIntegrationType->GetNumIntegrationPoints());
	if (mConstitutiveLaw!=0)
	{
		for (unsigned int count=0; count<mStaticData.size(); count++)
		{
			mStaticData[count] = rElement->AllocateStaticData(mConstitutiveLaw);
		}
	}
}
