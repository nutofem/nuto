// $Id: $
#include "nuto/mechanics/elements/ElementDataConstitutiveStaticData.h"
#include "nuto/mechanics/constitutive/ConstitutiveBase.h"
#include <assert.h>


NuTo::ElementDataConstitutiveStaticData::ElementDataConstitutiveStaticData(const NuTo::IntegrationTypeBase* rIntegrationType) : ElementDataBase(rIntegrationType)
{
	mConstitutiveLaw = 0;
}

NuTo::ElementDataConstitutiveStaticData::~ElementDataConstitutiveStaticData()
{
    for (int count=0; count< mStaticData.size(); count++)
    {
    	if (mStaticData[count]!=0)
    		delete mStaticData[count];
    }
}

void NuTo::ElementDataConstitutiveStaticData::SetConstitutiveLaw(const ElementBase* rElement, NuTo::ConstitutiveBase* rConstitutiveLaw)
{
	mConstitutiveLaw = rConstitutiveLaw;
	mStaticData.resize(mIntegrationType->GetNumIntegrationPoints());
	for (int count=0; count<mIntegrationType->GetNumIntegrationPoints(); count++)
	{
		try
		{
			mStaticData[count] = rElement->AllocateStaticData(mConstitutiveLaw);
		}
		catch (NuTo::MechanicsException &e)
		{
			e.AddMessage("[NuTo::ElementDataConstitutiveStaticData::SetConstitutiveLaw] Error allocating static data.");
			throw e;
		}
		catch (std::exception &e)
		{
			std::stringstream message;
			message << "[NuTo::ElementDataConstitutiveStaticData::SetConstitutiveLaw] std::exception " << e.what() << std::endl;
			throw MechanicsException(message.str());
		}
		catch (...)
		{
			throw MechanicsException("[NuTo::ElementDataConstitutiveStaticData::SetConstitutiveLaw] Error allocating static data.");
		}
	}
}

void NuTo::ElementDataConstitutiveStaticData::SetConstitutiveLaw(const ElementBase* rElement, int rIp, NuTo::ConstitutiveBase* rConstitutiveLaw)
{
	throw MechanicsException("[NuTo::ElementDataConstitutiveStaticData::SetConstitutiveLaw] Material can only be set for all integration points.");
}

NuTo::ConstitutiveStaticDataBase* NuTo::ElementDataConstitutiveStaticData::GetStaticData(int rIp)
{
	assert(rIp<(int)mStaticData.size());
	return mStaticData[rIp];
}

const NuTo::ConstitutiveStaticDataBase* NuTo::ElementDataConstitutiveStaticData::GetStaticData(int rIp)const
{
	assert(rIp<(int)mStaticData.size());
	return mStaticData[rIp];
}

NuTo::ConstitutiveBase* NuTo::ElementDataConstitutiveStaticData::GetConstitutiveLaw(int rIp)
{
	assert(rIp<(int)mStaticData.size());
	if (mConstitutiveLaw==0)
		throw MechanicsException("[NuTo::ElementDataConstitutiveStaticData::GetConstitutiveLaw] No material set for element.");
	return mConstitutiveLaw;
}

const NuTo::ConstitutiveBase* NuTo::ElementDataConstitutiveStaticData::GetConstitutiveLaw(int rIp)const
{
	assert(rIp<(int)mStaticData.size());
	if (mConstitutiveLaw==0)
		throw MechanicsException("[NuTo::ElementDataConstitutiveStaticData::GetConstitutiveLaw] No material set for element.");
	return mConstitutiveLaw;
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
